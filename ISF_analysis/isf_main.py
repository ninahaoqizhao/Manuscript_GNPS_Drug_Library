import os
import numpy as np
from pyteomics import mzml, mzxml
from numba import njit
import pandas as pd
import multiprocessing as mp
from tqdm import tqdm
import json


def main(drug_csv_path='data/analog_drug_pair_output.csv',
         msv_df_path='data/msv_results.tsv',
         mz_tol=0.01, rt_tol=10, extended_scans=10, out_dir='output',
         n_jobs=None):
    print('Creating data finder...')
    msv_dict = create_data_finder(msv_df_path)

    os.makedirs(out_dir, exist_ok=True)
    checkpoint_file = os.path.join(out_dir, 'checkpoint.json')

    # Load checkpoint data if it exists
    processed_files = {}
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as f:
            processed_files = json.load(f)
        print(f'Found checkpoint with {len(processed_files)} processed files')

    print('Reading and preprocessing data...')
    df = pd.read_csv(drug_csv_path)

    # Load existing output if it exists
    output_csv = os.path.join(out_dir, 'output.csv')
    if os.path.exists(output_csv):
        existing_df = pd.read_csv(output_csv)
        # Update the main dataframe with existing results
        df.update(existing_df)

    unique_files = df['analog_short_usi'].unique()
    # Filter out already processed files
    unprocessed_files = [f for f in unique_files if f not in processed_files]
    total_files = len(unprocessed_files)
    print(f'Processing {total_files} remaining files...')

    if total_files == 0:
        print('All files have been processed!')
        return

    process_args = []
    for file in unprocessed_files:
        sub_df = df[df['analog_short_usi'] == file].copy()
        sub_df = sub_df.drop_duplicates(subset=['drug_USI', 'analog_USI']).reset_index(drop=True)
        usi_pairs = [(drug, analog) for drug, analog in sub_df[['drug_USI', 'analog_USI']].values.tolist()]
        file_path = msv_dict.get(file.split('mzspec:')[1], None)

        if file_path:
            if not check_file_size(file_path):
                print(f'Skipping {file} due to large file size')
                continue
            process_args.append((file_path, usi_pairs, mz_tol, rt_tol, extended_scans, out_dir))

    if n_jobs is None:
        n_jobs = max(1, mp.cpu_count() - 1)
    print(f'Using {n_jobs} processes')

    files_processed = 0
    batch_results = {}

    with mp.Pool(processes=n_jobs) as pool:
        for i, (file, result) in enumerate(zip(unprocessed_files,
                                               tqdm(pool.imap_unordered(process_raw_data_file, process_args),
                                                    total=total_files,
                                                    desc='Processing files'))):

            if result is not None:
                # Update results in the dataframe
                for (drug_usi, analog_usi), ppc in result.items():
                    df.loc[(df['drug_USI'] == drug_usi) &
                           (df['analog_USI'] == analog_usi), 'ppc'] = ppc

                # Mark file as processed and store results
                processed_files[file] = True
                files_processed += 1

            # Save checkpoint every 100 files
            if files_processed > 0 and files_processed % 100 == 0:
                print(f'\nSaving checkpoint after {files_processed} files...')
                # Save checkpoint
                with open(checkpoint_file, 'w') as f:
                    json.dump(processed_files, f)
                # Save current results
                df.to_csv(output_csv, index=False)

    # Final save
    print('\nSaving final results...')
    with open(checkpoint_file, 'w') as f:
        json.dump(processed_files, f)
    df.to_csv(output_csv, index=False)
    print('Done!')


def check_file_size(file_path):
    file_size = os.path.getsize(file_path)
    # 1 Gb
    if file_size > 1e9:
        return False
    return True


def create_data_finder(msv_df_path):

    df = pd.read_csv(msv_df_path, sep='\t')
    df['msv_id'] = df['Filename'].apply(lambda x: x.split('/')[5])
    df['file'] = df['Filename'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    df['msv_file'] = df.apply(lambda x: f"{x['msv_id']}:{x['file']}", axis=1)

    # create a dictionary mapping from msv_file to Filename
    msv_dict = dict(zip(df['msv_file'], df['Filename']))

    return msv_dict


def process_raw_data_file(args):
    try:
        file_path, usi_pairs, mz_tol, rt_tol, extended_scans, out_dir = args
        file_format = os.path.splitext(file_path)[1].lower()

        usi_list = []
        for pair in usi_pairs:
            usi_list.append(pair[0])
            usi_list.append(pair[1])

        if file_format == '.mzml':
            ms1_scan_rts, targets = process_mzml(file_path, mz_tol, usi_list)
        elif file_format == '.mzxml':
            ms1_scan_rts, targets = process_mzxml(file_path, mz_tol, usi_list)
        else:
            raise ValueError(f"Unsupported file format: {file_format}")

        return calculate_ppc(usi_pairs, ms1_scan_rts, targets, rt_tol, extended_scans)
    except:
        return None


def process_mzml(file_path, mz_tol, usi_list):

    target_ms2_scans = [int(usi.split(':')[-1]) for usi in usi_list]

    # first generate target mzs, rts
    targets = []
    for spectrum in mzml.MzML(file_path):
        if spectrum['ms level'] == 2:  # MS2 scans only
            scan_number = spectrum['index'] + 1
            if scan_number not in target_ms2_scans:
                continue

            prec_mz = float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
            rt = spectrum['scanList']['scan'][0]['scan start time'] * 60  # in seconds

            targets.append({
                'usi': usi_list[target_ms2_scans.index(scan_number)],
                'ms2_scan': scan_number,
                'ms2_mz': prec_mz,
                'ms2_rt': rt,
                'feature_rt': None,
                'xic': []
            })

    # loop through again to get XICs
    ms1_scan_rts = []
    for spectrum in mzml.MzML(file_path):
        if spectrum['ms level'] == 1:

            rt = spectrum['scanList']['scan'][0]['scan start time'] * 60  # in seconds
            ms1_scan_rts.append(rt)

            mz_array = np.array(spectrum['m/z array'])
            intensity_array = np.array(spectrum['intensity array'])
            peaks = np.column_stack((mz_array, intensity_array))

            for target in targets:
                target['xic'].append(_get_peak_intensity(peaks, target['ms2_mz'], mz_tol))

    return ms1_scan_rts, targets


def process_mzxml(file_path, mz_tol, usi_list):

    target_ms2_scans = [int(usi.split(':')[-1]) for usi in usi_list]

    # first generate target mzs, rts
    targets = []
    for spectrum in mzxml.MzXML(file_path):
        if spectrum['msLevel'] == 2:  # MS2 scans only
            scan_number = spectrum['num']
            if scan_number not in target_ms2_scans:
                continue

            prec_mz = float(spectrum['precursorMz'][0]['precursorMz'])
            rt = spectrum['retentionTime']

            if isinstance(rt, str):
                rt = float(rt.replace("PT", "").replace("S", ""))
            else:
                rt = float(rt)
            rt = rt * 60  # in seconds

            targets.append({
                'usi': usi_list[target_ms2_scans.index(scan_number)],
                'ms2_scan': scan_number,
                'ms2_mz': prec_mz,
                'ms2_rt': rt,
                'feature_rt': None,
                'xic': []
            })

    # loop through again to get XICs
    ms1_scan_rts = []
    for spectrum in mzxml.MzXML(file_path):
        if spectrum['msLevel'] == 1:

            rt = spectrum['retentionTime']
            if isinstance(rt, str):
                rt = float(rt.replace("PT", "").replace("S", ""))
            else:
                rt = float(rt)
            rt = rt * 60

            ms1_scan_rts.append(rt)

            mz_array = np.array(spectrum['m/z array'])
            intensity_array = np.array(spectrum['intensity array'])
            peaks = np.column_stack((mz_array, intensity_array))

            for target in targets:
                target['xic'].append(_get_peak_intensity(peaks, target['ms2_mz'], mz_tol))

    return ms1_scan_rts, targets


def calculate_ppc(usi_pairs, ms1_scan_rts, targets, rt_tol=10, extended_scans=10):
    """
    Calculate Pearson correlation for drug-analog pairs
    """
    usi_pairs_ppc_dict = {}
    for drug_usi, analog_usi in usi_pairs:
        # Find corresponding targets
        drug = next((t for t in targets if t['usi'] == drug_usi), None)
        analog = next((t for t in targets if t['usi'] == analog_usi), None)

        if not drug or not analog:
            usi_pairs_ppc_dict[(drug_usi, analog_usi)] = 0.0
            continue

        # Get aligned XICs
        drug_xic, analog_xic = find_xics(ms1_scan_rts, drug, analog, rt_tol, extended_scans)

        # Calculate Pearson correlation
        drug_array = np.array(drug_xic)
        analog_array = np.array(analog_xic)

        if np.all(drug_array == 0) or np.all(analog_array == 0):
            correlation = 0.0
        else:
            correlation = np.corrcoef(drug_array, analog_array)[0, 1]
            if np.isnan(correlation):
                correlation = 0.0

        usi_pairs_ppc_dict[(drug_usi, analog_usi)] = correlation

    return usi_pairs_ppc_dict


def find_xics(ms1_scan_rts, drug, analog, rt_tol=10, extended_scans=10):
    """
    Extract aligned XICs for drug-analog pairs using drug's peak apex
    """
    ms1_scan_rts = np.array(ms1_scan_rts)
    xic_length = 2 * extended_scans + 1

    # Find peak apex using drug
    ms2_rt = drug['ms2_rt']
    min_rt = ms2_rt - rt_tol
    max_rt = ms2_rt + rt_tol

    idx = np.where((ms1_scan_rts >= min_rt) & (ms1_scan_rts <= max_rt))[0]
    if len(idx) == 0:
        return [0.0] * xic_length, [0.0] * xic_length

    intensities = np.array(drug['xic'])[idx]
    max_intensity_idx = idx[np.argmax(intensities)]

    # Extract aligned XICs
    start_idx = max_intensity_idx - extended_scans
    end_idx = max_intensity_idx + extended_scans + 1

    if start_idx < 0:
        padding = [0.0] * abs(start_idx)
        drug_xic = padding + drug['xic'][0:end_idx]
        analog_xic = padding + analog['xic'][0:end_idx]
    elif end_idx > len(ms1_scan_rts):
        padding = [0.0] * (end_idx - len(ms1_scan_rts))
        drug_xic = drug['xic'][start_idx:] + padding
        analog_xic = analog['xic'][start_idx:] + padding
    else:
        drug_xic = drug['xic'][start_idx:end_idx]
        analog_xic = analog['xic'][start_idx:end_idx]

    return drug_xic[:xic_length], analog_xic[:xic_length]


@njit
def _get_peak_intensity(peaks, mz, tol):

    idx = np.abs(peaks[:, 0] - mz) < tol
    if not np.any(idx):
        return 0.0

    return np.max(peaks[idx, 1])


if __name__ == '__main__':

    # df = pd.read_csv('data/analog_drug_pair_output.csv')
    # print(df['analog_short_usi'].nunique())

    import argparse
    parser = argparse.ArgumentParser(description='Process drug-analog pairs')
    parser.add_argument('--test', action='store_true', help='Run test')
    parser.add_argument('--n_jobs', type=int, default=None, help='Number of processes to use')

    args = parser.parse_args()

    if args.test:
        main(drug_csv_path='data/test.csv',
             msv_df_path='data/msv_results.tsv',
             mz_tol=0.005, rt_tol=10, extended_scans=5, out_dir='output', n_jobs=args.n_jobs)
    else:
        main(drug_csv_path='data/full_drug_analog_pairs_filtered_masst_results.csv',
             msv_df_path='data/msv_results.tsv',
             mz_tol=0.005, rt_tol=10, extended_scans=5, out_dir='output', n_jobs=args.n_jobs)

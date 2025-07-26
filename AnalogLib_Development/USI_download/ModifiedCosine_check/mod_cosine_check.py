import numpy as np
from matchms import Spectrum
from matchms.similarity import ModifiedCosine
import pandas as pd
from tqdm import tqdm


def read_mgf_to_dict(library_mgf):
    """
    id: spec
    """

    allowed_keys = ['TITLE', 'PEPMASS', 'SPECTRUMID']

    with open(library_mgf, 'r') as file:
        all_spec = {}
        for line in file:
            # empty line
            _line = line.strip()
            if not _line:
                continue
            elif line.startswith('BEGIN IONS'):
                spectrum = {}
                # initialize spectrum
                mz_list = []
                intensity_list = []
            elif line.startswith('END IONS'):
                if len(mz_list) == 0:
                    continue
                spectrum['mz_ls'] = mz_list
                spectrum['intensity_ls'] = intensity_list

                all_spec[spectrum['id']] = spectrum

                continue
            else:
                # if line contains '=', it is a key-value pair
                if '=' in _line:
                    # split by first '='
                    key, value = _line.split('=', 1)

                    if key not in allowed_keys:
                        continue

                    if key == 'PEPMASS':
                        spectrum['prec_mz'] = value
                    else:
                        if 'CCMSLIB' in value or 'mzspec:' in value:
                            spectrum['id'] = value
                else:
                    # if no '=', it is a spectrum pair
                    this_mz, this_int = _line.split()
                    try:
                        mz_list.append(float(this_mz))
                        intensity_list.append(float(this_int))
                    except:
                        continue

    return all_spec


def main(input_csv, input_mgf, mz_tol=0.05):
    mod_cos_eng = ModifiedCosine(tolerance=mz_tol)

    df = pd.read_csv(input_csv)

    # Read the MGF file into a DataFrame
    mgf_dict = read_mgf_to_dict(input_mgf)

    df['score'] = None
    df['matched_peaks'] = None

    for i, row in tqdm(df.iterrows(), total=len(df)):
        analog_spec = mgf_dict.get(row['analog'], None)
        ref_spec = mgf_dict.get(row['drug'], None)

        if analog_spec is None or ref_spec is None:
            continue

        try:
            # Convert to matchms Spectrum objects
            analog_spec_obj = Spectrum(mz=np.array(analog_spec['mz_ls']),
                                       intensities=np.array(analog_spec['intensity_ls']),
                                       metadata={"precursor_mz": analog_spec['prec_mz']})
            ref_spec_obj = Spectrum(mz=np.array(ref_spec['mz_ls']),
                                    intensities=np.array(ref_spec['intensity_ls']),
                                    metadata={"precursor_mz": ref_spec['prec_mz']})

            # Calculate the modified cosine similarity
            score = mod_cos_eng.pair(analog_spec_obj, ref_spec_obj)

            df.at[i, 'score'] = score['score']
            df.at[i, 'matched_peaks'] = score['matches']
        except:
            print('Error:', row['analog'], row['drug'])
            continue

    df.to_csv('output.csv', index=False)


if __name__ == '__main__':

    main('drug_analog_usi_pairs.csv', 'input_all_spectra_drug_analog_usi.mgf', mz_tol=0.05)

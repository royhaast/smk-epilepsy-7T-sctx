function b1correct(out_folder)
    warning('off','all')

    addpath(genpath('./code_cemerem'));
    addpath(genpath('../presurfer/spm12'));
    
    MP2RAGE_Processing_Magnetom_BIDS(out_folder)

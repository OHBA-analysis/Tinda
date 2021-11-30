conf = mvarconfig('mvarconfig.ini');
successfulruns = false(length(conf.PARTICIPANTS),1);

for i=1:length(conf.PARTICIPANTS)
    try
        mvar_02_hcp_aal_beamform_CHedit(conf.PARTICIPANTS{i});
        successfulruns(i) = true;
    catch ME
        fprintf('Failed')
        successfulruns(i) = false;
        errormessages{i} = ME.message;
    end
end 


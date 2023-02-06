function [info, varnames, age, sex, RT] = getparticipantinfo(whichstudy)

config = getStudyDetails(whichstudy);

if whichstudy==1
  [info,varnames] = MEGUK_getparticipantinfo();
  age = info(:,11);
  sex = info(:,12); % 1 denotes male, 2 denotes female
  RT = info(:,3);
elseif whichstudy==3
  [info,varnames] = HCP_getparticipantinfo(config);
  age = info(:,4);
  sex = info(:,3); % 1 denotes female
  RT = info(:,end);
elseif whichstudy==4
  info = camcan_getparticipantinfo(config);
  age = info(:,1);
  sex = info(:,3); % 1 denotes male, 2 female
  RT = info(:,5);
end
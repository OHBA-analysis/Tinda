function color_scheme = colorscheme(whichstudy)

if whichstudy==1 || whichstudy==2
  % Colourscheme based on the ColorBrewer Set 1
  
  % http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=8
  set1 = { [228,26,28], [55,126,184], [77,175,74], [152,78,163],...
    [255,127,0], [247,129,191], [255,255,51],[166,86,40]};
  
  %K=12 option:
  set1 = {[141,211,199],[188,128,189],[190,186,218],[251,128,114], ...
    [128,177,211], [253,180,98],[179,222,105],[252,205,229], ...
    [217,217,217],[255,255,179],[204,235,197],[255,237,111]};
  
  
  % normalise and convert back to cell array
  set1 = cellfun(@(x) (0.9*x)./255, set1,'UniformOutput', false );
  
  % actually better just to use inbuilt matlab colors:
  colors_alt=jet(12);
  rng(10);
  colors_alt = colors_alt(randperm(12),:);
  for k=1:12
    set1{k}=colors_alt(k,:);
  end
  set1{8} = set1{8} + 0.5*[0,1,1];
  set1{7} = [141,211,199]./255;
  set1{11} = [190,186,218]./255;
  color_scheme = set1;
elseif whichstudy==3 || whichstudy==5
  tmp = circshift([0.400000000000000,0.760784313725490,0.647058823529412;0.988235294117647,0.552941176470588,0.384313725490196;0.552941176470588,0.627450980392157,0.796078431372549;0.905882352941177,0.541176470588235,0.764705882352941;0.650980392156863,0.847058823529412,0.329411764705882;1,0.850980392156863,0.184313725490196;0.898039215686275,0.768627450980392,0.580392156862745;0.701960784313725,0.701960784313725,0.701960784313725;0.400000000000000,0.260784313725490,0.647058823529412;0.738235294117647,0.352941176470588,0.384313725490196;0.552941176470588,0.827450980392157,0.946078431372549;0.605882352941177,0.541176470588235,0.164705882352941],-2,1);
  for k=1:12
    color_scheme{k} = tmp(k,:);
  end
elseif whichstudy==4 || whichstudy==6 || whichstudy==8
  tmp = brewermap(12, 'Set3');
  for k=1:12
    color_scheme{k} = tmp(k,:);
  end
elseif whichstudy==7
      tmp = colormap(hsv);
  for k=1:12
    color_scheme{k} = tmp(1+(k-1)*floor(size(tmp,1)/12),:);
  end
end
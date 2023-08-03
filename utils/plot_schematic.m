
color_scheme = set1_cols();
T1 = [0, 4, 7, 3, 2, 10, 2, 9, 8, 4, 5, 2];
T2 = [0, 3, 4, 0, 10, 6, 5, 4, 8, 5, 2, 2];

%% Bar graph
figure; subplot(1,2,1), hold on
for ii=1:12
  h=barh(ii, T1(ii));
  set(h,'FaceColor', color_scheme{ii})
end
set(gca, 'Ydir', 'reverse')
set(gca, 'YTick', 1:12)
subplot(1,2,2), hold on
for ii=1:12
  h=barh(ii, T2(ii));
  set(h,'FaceColor', color_scheme{ii})
end
set(gca, 'Ydir', 'reverse')
set(gca, 'YTick', 1:12)

%% 
q = squeeze(FO(1,:,:,:));
clear d

cb = [256,193,1; 201, 204, 231]/256;
ii=2;
figure;
nrand=100;
for ii=2:12
  subplot(12,1,ii),
  if ii==5 || ii==7
    d{1} = normrnd(0,1,[1,nrand]);
    d{2} = normrnd(0.2*ii+rand(1),1,[1,nrand]);
  elseif ii==2 || ii==3 || ii==4 || ii==6
    d{1} = normrnd(0.2*ii+rand(1),1,[1,nrand]);
    d{2} = normrnd(0,1,[1,nrand]);
  else
    d{1} = normrnd(0,1,[1,nrand]);
    d{2} = normrnd(0,1,[1,nrand]);
  end
% d{1} = %squeeze(FO(1,ii,1,:));
% d{2} = %squeeze(FO(1,ii,2,:));

h1 = raincloud_plot(d{1}, 'box_on', 0, 'color', cb(1,:), 'alpha', 0.5);
h2 = raincloud_plot(d{2}, 'box_on', 0, 'color', cb(2,:), 'alpha', 0.5);
ylim([0 .7])
% view([90 -90])
axis off
end

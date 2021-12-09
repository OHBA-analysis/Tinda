function fig1 = plotCyclicalTimeFreqPattern(manualorder,p,a,b,f)
% state sequence is the sequence of states, starting at the noon position
% and proceeding clockwise

set(gcf,'Position',[1 54 1440 744]);
applythresh = true; % threshold spatial maps

[K,nch,num_nodes] = size(a);
% preliminary: setup rotational geometry of states:
disttoplot_manual = zeros(12,2);
for i=1:12
    %temp = exp(sqrt(-1)*((-i+1)/12*2*pi+pi/2));
    temp = exp(sqrt(-1)*(i+2)/12*2*pi);
    disttoplot_manual(manualorder(i),:) = [real(temp),imag(temp)];
end
color_scheme = set1_cols;

subplot(1,3,1);

for ik=1:K
    scatter1 = scatter(disttoplot_manual(ik,1),disttoplot_manual(ik,2),400,...
        'MarkerFaceColor',color_scheme{ik},'MarkerEdgeColor',color_scheme{ik}); 
    hold on
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatter1.MarkerFaceAlpha = 1;%.75;
    
    if ik<10
        text(disttoplot_manual(ik,1)-0.01,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    else
        text(disttoplot_manual(ik,1)-0.07,disttoplot_manual(ik,2),int2str(ik),'FontSize',12,'FontWeight','bold');hold on;
    end
end
axis square
axis off
vidres = 1000;
freq_time_map = cell(num_nodes,1);
weights = zeros(vidres,K);
for t=1:vidres
    temp = exp(sqrt(-1)*((-t-1)/vidres*2*pi + pi/2));
    xpdf = mvnpdf(disttoplot_manual,[real(temp),imag(temp)],0.05*eye(2));
    weights(t,:) = xpdf./sum(xpdf);
    for inode = 1:num_nodes
        freq_time_map{inode}(:,t) = weights(t,:)*a(:,:,inode);
    end
end

% and plot:
%nodelocs = setdiff(1:9,1:3:9);
nodelocs = 3:3:(6*6);
freq_labels = [10:10:max(f)];
for i=1:length(freq_labels)
    freq_locs(i) = find(f>freq_labels(i),1);
end
freq_locs = fliplr(length(f)-freq_locs+1);
for inode = 1:6
    subplot(6,3,nodelocs(inode));
    %make time-frequency plot:
    freq_time_map{inode} = demean(freq_time_map{inode},2);
    imagesc(flipud(freq_time_map{inode}));
    set(gca,'YTick',freq_locs);
    set(gca,'YTickLabel',fliplr(freq_labels));
    set(gca,'XTick',[]);
    plot4paper('Cycle time','Frequency');
end
colormap('hot');
cm = colormap;
colormap(cm(1:50,:));
% plot 2 way brain plots:
surface_inflation = 2;
interptype = 'trilinear';
for inode=1:6
    data = b(inode,:);
    
    if applythresh
        data(data<prctile(data,80))=NaN;
    end
    niifile = p.savenii(data);
    output_right    = [niifile '_right.func.gii'];
    output_left     = [niifile '_left.func.gii'];
    cl = onCleanup(@()  cellfun(@delete,{niifile,output_left,output_right})); % Enable deleting temp files even if debugging

    surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii');
    surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii');

    switch surface_inflation
        case 0
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii');
        case 1
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.inflated.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.inflated.32k_fs_LR.surf.gii');
        case 2
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.very_inflated.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.very_inflated.32k_fs_LR.surf.gii');
    end
    runcmd('wb_command -volume-to-surface-mapping %s %s %s -%s',niifile,surf_right,output_right,interptype)
     runcmd('wb_command -volume-to-surface-mapping %s %s %s -%s',niifile,surf_left,output_left,interptype)

    sl = gifti(display_surf_left);
    vl = gifti(output_left);
    sr = gifti(display_surf_right);
    vr = gifti(output_right);
    %hfig = figure('Position',[547 100 577 453]);
    hfig = gcf;
    set(hfig,'Color','White');
    if isempty([])
        clims = [min([min(vl.cdata) min(vr.cdata)]) max([max(vl.cdata) max(vr.cdata)])];
    end
    vl.cdata(vl.cdata==0)=NaN;
    vr.cdata(vr.cdata==0)=NaN;
    % make colormap
    if min([min(vl.cdata) min(vr.cdata)]) < 0 || clims(1) < 0
        cm = cat(2,linspace(.5, 1 ,63)',linspace(0, 1 ,63)',linspace(0, 0 ,63)');
        cm2 = cat(2,linspace(0, 0 ,63)',linspace(1, 0 ,63)',linspace(1, .5 ,63)');
        %cm = cat(1,cm2,[.6 .6 .6],cm);
        cm = cat(1,flipud(cm2),flipud(cm));
    else
        cm = cat(2,linspace(.5, 1 ,63)',linspace(0, 1 ,63)',linspace(0, 0 ,63)');
        c=flipud(cm);
        cm = cat(1,[.6 .6 .6],cm);
    end
    
    for ii = 1:4
%         if ii==3
%             ax(ii) = axes('Position',[.05+((ii-1)*.165) .1 .175 .8]);
%         else
%             ax(ii) = axes('Position',[.05+((ii-1)*.165) .3 .175 .5]);
%         end
        %ax(ii) = subplot(2,2,ii);
    end
%     ax(1) = axes('Position',[0 0.5 0.5 0.5]);
%     ax(2) = axes('Position',[0.55 0.5 0.5 0.5]);
%     ax(3) = axes('Position',[0.09 0.01 0.5 0.5]);
%     ax(4) = axes('Position',[0.475 0.01 0.5 0.5]);
    
    ax(1) = subplot(6,6,3+(inode-1)*6);
    ax(2) = subplot(6,6,4+(inode-1)*6);
    %ax(3) = subplot(6,6,6+3);
    %ax(4) = subplot(6,6,6+4);
    ax(1).Position(3:4) = 1.2*ax(1).Position(3:4);
    ax(2).Position(3:4) = 1.2*ax(2).Position(3:4);
    ax(1).Position(4) = 1.2*ax(1).Position(4);
    ax(2).Position(4) = 1.2*ax(2).Position(4);
    
    % lateral
%     axes(ax(3));
%     s(1) = patch('Faces',sl.faces,'vertices',sl.vertices,'CData',[]);
%     hold on
%     sg(1) = patch('Faces',sl.faces,'vertices',sl.vertices);
%     
%     axes(ax(4));
%     s(2) = patch('Faces',sr.faces,'vertices',sr.vertices,'CData',[]);
%     hold on
%     sg(2) = patch('Faces',sr.faces,'vertices',sr.vertices);    
%     
%     set(s(1),'FaceVertexCData',vl.cdata)
%     set(s(2),'FaceVertexCData',vr.cdata)
% 
%     set(sg(1),'FaceVertexCData',0.4*ones(size(vl.cdata,1),3));
%     set(sg(2),'FaceVertexCData',0.4*ones(size(vr.cdata,1),3));
%     set(sg(1),'FaceVertexAlphaData',+~isfinite(vl.cdata),'FaceAlpha','interp','AlphaDataMapping','none');
%     set(sg(2),'FaceVertexAlphaData',+~isfinite(vr.cdata),'FaceAlpha','interp','AlphaDataMapping','none');    
    
    view(ax(1),[270 0])
    view(ax(2),[-270 0])
    clear s sg
    
    % medial:
    axes(ax(1));
    s(1) = patch('Faces',sl.faces,'vertices',sl.vertices,'CData',[]);
    hold on
    sg(1) = patch('Faces',sl.faces,'vertices',sl.vertices);
    
    axes(ax(2));
    s(2) = patch('Faces',sr.faces,'vertices',sr.vertices,'CData',[]);
    hold on
    sg(2) = patch('Faces',sr.faces,'vertices',sr.vertices);    
    
    set(s(1),'FaceVertexCData',vl.cdata)
    set(s(2),'FaceVertexCData',vr.cdata)

    set(sg(1),'FaceVertexCData',0.4*ones(size(vl.cdata,1),3));
    set(sg(2),'FaceVertexCData',0.4*ones(size(vr.cdata,1),3));
    set(sg(1),'FaceVertexAlphaData',+~isfinite(vl.cdata),'FaceAlpha','interp','AlphaDataMapping','none');
    set(sg(2),'FaceVertexAlphaData',+~isfinite(vr.cdata),'FaceAlpha','interp','AlphaDataMapping','none');    
    
%     view(ax(3),[90 0])
%     view(ax(4),[-90 0])
    clear s sg
    
    
    arrayfun(@(x) shading(x,'interp'),ax);
    arrayfun(@(x) axis(x,'off'), ax);

    
    for ii = 1:2%4
        axes(ax(ii))
        %colormap(cm)
        %zoom(1.1)
        if ~isempty(clims)
            caxis(clims);
        else
            mn = min([min(vl.cdata) min(vr.cdata)]);
            mx = max([max(vl.cdata) max(vr.cdata)]);
            if mn < 0
                m = max([ abs(mn) abs(mx) ]);
                caxis([-m m]);
            else
                caxis([0 mx])
            end
        end
        material( [0.3, 0.8, 0.2] );
        camlight('left')
        camlight('right')
        camlight('headlight')
        if ii == 5
            c = colorbar('EastOutside');
            c.Position(1) = .9;
            c.Position(2) = .25;
            c.Position(3) = .02;
            c.FontSize = 18;
        end
    end
end
end

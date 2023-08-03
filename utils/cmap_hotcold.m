function hotcold = cmap_hotcold

cm = colormap(hot);close
cm2 = cm;
cm2(:,1) = cm(:,3);
cm2(:,3) = cm(:,1);
hotcold=[flipud(cm2);cm];
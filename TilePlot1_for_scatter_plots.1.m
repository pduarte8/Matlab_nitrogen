prompt = {'Number of lines:','Number of columns:'};
dlg_title = 'Tilelayout dimensions';
num_lines= 1;
def = {'2','2'}
answer  = inputdlg(prompt,dlg_title,num_lines,def);
n = str2double(answer(1));
m = str2double(answer(2));

Fig1 = figure;

t = tiledlayout(n,m,'TileSpacing','compact','Padding','compact');

for i = 1:n*m
   [filename, filepath] = uigetfile('*.fig', 'Open File')
   fullname = fullfile(filepath, filename);
   existingFig = openfig(fullname);
   aax = gca
   %Ticks = aax.XTick
   %Labels = aax.XTickLabel
   figure(Fig1);
   axNew = nexttile(i);
   copyobj(allchild(get(existingFig,'CurrentAxes')),axNew);

% Copy axis properties
   axNew.XTick = aax.XTick;
   axNew.XTickLabel = aax.XTickLabel
   axNew.XLim = aax.XLim;
   axNew.YLim = aax.YLim;
   axNew.XLabel.String = aax.XLabel.String;
   axNew.YLabel.String = aax.YLabel.String;
   axNew.Title.String  = aax.Title.String;
   axNew.Box = aax.Box;
   axNew.GridLineStyle = aax.GridLineStyle;
   grid(axNew,'on')
    
   % --- Preserve legend ---
   lgdOld = findobj(aax,'Type','legend');
   if ~isempty(lgdOld)
        legend(axNew, ...
            lgdOld.String, ...
            'Location', lgdOld.Location);
   end

   %set(gca,'fontsize',10)
   %prompt = {'Y axis lower limit','Y axis upper limit'};
   %dlg_title = 'Axis limits';
   %num_lines= 1;
   %def = {'0','10'}
   %answer  = inputdlg(prompt,dlg_title,num_lines,def);
   %ylim([str2double(answer(1)) str2double(answer(2))]);

  
   close(existingFig);
end  
t.TileSpacing="Compact"
t.Padding="Compact"
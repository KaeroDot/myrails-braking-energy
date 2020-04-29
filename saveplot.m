%% --- Function saveplot --------------------
% saves plot into gnuplot script if octave is used, or into png if matlab
function saveplot(filename, dirpath)
        if (exist ("OCTAVE_VERSION", "builtin") > 0)
                printplt(fullfile(dirpath, filename));
        else
                saveas(gcf,fullfile(dirpath, [filename '.fig']))
        end
end

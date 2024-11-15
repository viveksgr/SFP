function SFP_resetfigs()
% Set default font to Arial for all text in figures
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultUicontrolFontName', 'Arial');   % Affects UI controls like buttons
set(0, 'DefaultUipanelFontName', 'Arial');     % Affects UIPanels
set(0, 'DefaultUitableFontName', 'Arial');     % Affects UITables
end

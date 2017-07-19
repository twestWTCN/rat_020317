function mean_diffstat_tabsave(stattab,feat,R)
delete([R.analysispath R.pipestamp '\results\datatables\csv\mean_diffstat_' feat '.csv'])
% This can easily be adapted to export any table!!
sourcetab = vertcat(stattab{:});
if ~exist([R.analysispath R.pipestamp '\results\datatables\matlab\'], 'dir')
    mkdir([R.analysispath R.pipestamp '\results\datatables\matlab\']);
end
save([R.analysispath R.pipestamp '\results\datatables\matlab\mean_diffstat_' feat '.mat'],'sourcetab')

% Save text version
if ~exist([R.analysispath R.pipestamp '\results\datatables\csv\'], 'dir')
    mkdir([R.analysispath R.pipestamp '\results\datatables\csv\']);
end
cd([R.analysispath R.pipestamp '\results\datatables\csv'])
fid = fopen(['mean_diffstat_' feat '.csv'],'a');
dlmwrite(['mean_diffstat_' feat '.csv'],sourcetab,'-append','precision','%.6f','delimiter',',');
fid = fopen(['mean_diffstat_' feat '.csv'],'r'); fclose('all');

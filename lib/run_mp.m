function task = run_mp(tag, m, mpopt, fname, solvedcase)

if nargin < 5
    solvedcase =  '';
    if nargin < 4
        fname = '';
        if nargin < 3
            mpopt = mpoption;
        end
    end
end

%% create task object
switch upper(tag)
    case 'PF'
        mp_task_class = @mp_task_pf;
    case 'CPF'
        mp_task_class = @mp_task_cpf;
    case 'OPF'
        mp_task_class = @mp_task_opf;
end
task = mp_task_class();

%% run task
task.run(m, mpopt);

%% pretty-print results to console & possibly to file
task.print_soln(fname, mpopt);

%% save solved case
if ~isempty(solvedcase) && task.success
    task.save_soln(solvedcase);
end

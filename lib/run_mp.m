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
switch upper(tag)
    case 'PF'
        mp_task_class = @mp_task_pf;
    case 'CPF'
        mp_task_class = @mp_task_cpf;
    case 'OPF'
        mp_task_class = @mp_task_opf;
end

task = mp_task_class();

dm = task.create_data_model(m, mpopt);
nm = task.create_network_model(dm, mpopt);
mm = task.create_math_model(nm, dm, mpopt);

success = task.run(mm, nm, dm, mpopt);

if success
    %% update network model with math model solution
    task.mm2nm(mm, nm);

    %% update data model with network model solution
    task.nm2dm(nm, dm, mpopt);

    %% pretty-print results to console & possibly to file
    task.print_soln(fname);

    %% save solved case
    if ~isempty(solvedcase)
        task.save_soln(solvedcase);
    end
end
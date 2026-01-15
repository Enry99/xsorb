'''
For each calc_type, print the status of all calculations in the database-

'''

from xsorb.adsorptiondata.adsorptioncalculation import ALLOWED_STATUSES
from xsorb.io.database import Database
from xsorb.io.jobs import get_running_and_queued_job_ids


def print_status(calc_type: str) -> None:
    '''
    For each calc_type, print the status of all calculations in the database,
    in the format:

    Screening:
        completed: [1,2,5...]
        incomplete: [3,4(r),6...]
        scf_nonconverged: [6,7(q)...]
    Relax:
        completed: [8,9...]
        incomplete: [10(r),11...]
        scf_nonconverged: [6,7(q)...]
    Mlopt:
        completed: [12...]
        incomplete: [13(q),14...]

    where (r) indicates a running job, and (q) a queued job.
    '''

    print("Current status of Xsorb calculations:")

    calc_types = []
    if calc_type == 'all':
        calc_types = ['screening', 'relax', 'mlopt']
    else:
        calc_types = [calc_type]

    running_job_ids, queued_job_ids = get_running_and_queued_job_ids()

    for ctype in calc_types:
        print(f"{ctype.capitalize()}:")

        #update database
        Database.update_calc_db(calc_type=ctype)

        for status in ALLOWED_STATUSES:
            rows = Database.get_calculations(calc_type=ctype,
                                             selection=f'status={status}',
                                             update=False,
                                             verbose=False)
            ids = []
            for row in rows:
                jid = row.job_id
                if jid in running_job_ids:
                    ids.append(f"{row.calc_id}(r)")
                elif jid in queued_job_ids:
                    ids.append(f"{row.calc_id}(q)")
                else:
                    ids.append(f"{row.calc_id}")
            print(f"    {status}: {', '.join(ids)}")
        #print("")
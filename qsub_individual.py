import re

REGEX = re.compile('ds_dataset(?P<dataset_id>\d)_\d+-(?P<weight>\d\.?\d?)-(?P<support>\w+)-(?P<gridifier>\w+)-(?P<pipeline>\w+)-(?P<run>\d)')

ds = [
    "ds_dataset3_80-0-path-hagrid-heuristic-4",
    "ds_dataset3_80-0-path-hagrid-heuristic-3",
    "ds_dataset3_80-0-path-hagrid-heuristic-2",
    "ds_dataset3_80-0-path-hagrid-heuristic-5",
    "ds_dataset3_80-0-path-hagrid-heuristic-1",
    "ds_dataset5_120-0.5-path-dgrid-heuristic-3",
    "ds_dataset5_120-0.5-path-dgrid-heuristic-4",
    "ds_dataset5_120-0.5-path-dgrid-heuristic-5",
    "ds_dataset5_120-0.5-path-dgrid-heuristic-2",
    "ds_dataset5_120-0.5-path-dgrid-heuristic-1",
    "ds_dataset7_160-0.5-path-hagrid-heuristic-5",
    "ds_dataset7_160-0.5-path-hagrid-heuristic-2",
    "ds_dataset7_160-0.5-path-hagrid-heuristic-3",
    "ds_dataset7_160-0.5-path-hagrid-heuristic-4",
    "ds_dataset7_160-0.5-path-hagrid-heuristic-1",
    "ds_dataset9_200-0.5-path-dgrid-heuristic-1",
    "ds_dataset9_200-0.5-path-dgrid-heuristic-4",
    "ds_dataset9_200-0.5-path-dgrid-heuristic-3",
    "ds_dataset9_200-0.5-path-dgrid-heuristic-2",
    "ds_dataset9_200-0.5-path-dgrid-heuristic-5",
]

def to_qsub(d):
    num_glyphs = (d['dataset_id'] + 1) * 20
    jobname = f"esvis-ds_dataset{d['dataset_id']}_{num_glyphs}-{d['weight']}-{d['support']}-{d['gridifier']}-{d['pipeline']}-{d['run']}"
    mem = '32G'
    dataset = f"ds_dataset{d['dataset_id']}"
    base = '/home1/npiccolotto/ensemble-sets/results'
    hlr = 18000
    return f"qsub -N {jobname} -l bc5 -l mem_free={mem} -l h_vmem={mem} -l h_rt={hlr} -e {base}/logs/ -o {base}/logs/ -r y run.sh {base}/{jobname} {d['support']} {d['gridifier']} {d['pipeline']} {dataset} {weight}"

for d in ds:
    match = REGEX.match(d)
    dataset_id = int(match.group('dataset_id'))
    weight = match.group('weight')
    support = match.group('support')
    gridifier = match.group('gridifier')
    pipeline = match.group('pipeline')
    run = int(match.group('run'))

    print(to_qsub({
        'dataset_id': dataset_id,
        'weight': weight,
        'support': support,
        'gridifier': gridifier,
        'pipeline': pipeline,
        'run': run
    }))


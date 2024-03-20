import dimod
import csv
from tqdm import tqdm
from pathlib import Path
from dwave.system import LeapHybridCQMSampler
from argparse import ArgumentParser


class Instance:
    def __init__(self, resource_caps, tasks=None, name="", max_length=None):
        if tasks is None:
            self.tasks: list[Task] = []
        else:
            self.tasks: list[Task] = tasks
        self.resource_caps = resource_caps
        self.name = name

        self.max_length = max_length
        # if max_length is None:
        #     self.max_length = sum([t.length for t in self.tasks])
        # else:
        #     self.max_length = max_length

    def calculate_bounds(self):
        if self.max_length is None:
            self.max_length = sum([t.length for t in self.tasks])
        # calculate lower bounds
        for task in self.tasks:
            for i in task.successors:
                self.tasks[i - 1].lower_bound = max(
                    self.tasks[i - 1].lower_bound, task.lower_bound + task.length
                )
        # calculate upper bounds
        for task in reversed(self.tasks):
            task.upper_bound = self.max_length - task.length
            for i in task.successors:
                task.upper_bound = min(
                    task.upper_bound, self.tasks[i - 1].upper_bound - task.length
                )

        # # calculate lower bounds
        # for task in self.tasks:
        #     for i in task.successors:
        #         self.tasks[i - 1].lower_bound = max(
        #             self.tasks[i - 1].lower_bound,
        #             task.lower_bound + task.length
        #         )
        # # calculate upper bounds
        # for task in reversed(self.tasks):
        #     task.upper_bound = self.max_length - task.length
        #     for i in task.successors:
        #         task.upper_bound = min(
        #             task.upper_bound,
        #             self.tasks[i - 1].upper_bound - task.length
        #         )

    def __repr__(self) -> str:
        res = f"PSP Instance: {self.name}\nnumber of tasks: {len(self.tasks)}\nnumber of resources: {len(self.resource_caps)}, max_length: {self.max_length}"
        return res

    def get_cqm(self):
        cqm = dimod.ConstrainedQuadraticModel()
        v = []
        task_vars = {}
        for task in self.tasks:
            v.append(
                dimod.Integer(
                    "t" + str(task.number),
                    lower_bound=task.lower_bound,
                    upper_bound=task.upper_bound,
                )
            )
            for t in range(task.lower_bound, task.upper_bound + 1):
                task_vars[(task.number, t)] = dimod.Binary(
                    f"task_{task.number}_time_{t}"
                )

        # setting up binary variables to be one-hot of integer
        for task in self.tasks:
            cqm.add_constraint(
                dimod.quicksum(
                    [
                        task_vars[(task.number, t)] * t
                        for t in range(task.lower_bound, task.upper_bound + 1)
                    ]
                )
                - v[task.number - 1]
                == 0,
                label=f"task {task.number} binary value",
            )
            cqm.add_constraint(
                dimod.quicksum(
                    [
                        task_vars[(task.number, t)]
                        for t in range(task.lower_bound, task.upper_bound + 1)
                    ]
                )
                == 1,
                label=f"task {task.number} one-hot",
            )

        # precedence constraint
        for task in self.tasks:
            t = task.number
            for s in task.successors:
                cqm.add_constraint(
                    v[s - 1] - v[t - 1] >= task.length, label=f"{t} preceeds {s}"
                )

        # for time in self.max_length:
        # cqm.add_constraint(
        #     3 <= dimod.quicksum([int(var <= 3) for var in v]), label='no more than 3 res between 1 and 3'
        # )

        # resource constraint
        for ir, r in enumerate(self.resource_caps):
            for timestamp in range(self.max_length + 1):
                sum = 0
                for task in self.tasks:
                    if not (task.lower_bound <= timestamp <= task.upper_bound):
                        continue
                    for t in range(
                        max(task.lower_bound, timestamp - task.length + 1),
                        min(task.upper_bound + 1, timestamp + 1),
                    ):
                        sum += task_vars[(task.number, t)] * task.resource_costs[ir]
                if type(sum) is int:
                    continue
                cqm.add_constraint(
                    sum <= r,
                    label=f"resource {ir} constraint for timestamp {timestamp}",
                )

        cqm.set_objective(v[-1])
        return cqm


class Task:
    def __init__(self, number, costs, length, successors):
        self.number = number
        self.resource_costs = costs
        self.length = length
        self.successors = successors
        self.lower_bound = 0
        self.upper_bound = 9999999

    def __repr__(self) -> str:
        divider_down = "_" * 30 + "\n"
        res = f"Task number {self.number}, length: {self.length}\nresource costs: {self.resource_costs}\nsuccessors: {self.successors}\nlower bound: {self.lower_bound}, upper bound: {self.upper_bound}"
        divider_up = "\n" + "-" * 30
        return divider_down + res + divider_up


def read_rcp_instance(filepath):
    with open(filepath, "r") as file:
        first_line = ""
        while not first_line.strip():
            first_line = next(file)  # no need to read number of resources and tasks
        resource_caps = list(map(int, next(file).strip().split()))
        num_resources = len(resource_caps)
        instance = Instance(resource_caps, name=Path(filepath).name)
        for line in file:
            line = list(map(int, line.strip().split()))
            if not line:
                continue
            length = line[0]
            resource_costs = line[1 : num_resources + 1]
            successors = line[num_resources + 2 :]
            instance.tasks.append(
                Task(len(instance.tasks) + 1, resource_costs, length, successors)
            )
        return instance


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--instance", type=str, required=True)
    parser.add_argument("--time_limit", type=int, default=None)
    parser.add_argument("--stats_file", type=str, default="stats.csv")
    args = parser.parse_args()
    return args


instances = [
    "./data/DC2/DC2/npv25/rcpspdc18.rcp",
    "./data/DC2/DC2/npv25/rcpspdc101.rcp",
    "./data/DC2/DC2/npv50/rcpspdc247.rcp",
    "./data/DC2/DC2/npv50/rcpspdc341.rcp",
    "./data/DC2/DC2/npv75/rcpspdc392.rcp",
    "./data/DC2/DC2/npv75/rcpspdc436.rcp",
    "./data/DC2/DC2/npv75/rcpspdc505.rcp",
    "./data/DC2/DC2/npv75/rcpspdc536.rcp",
    "./data/DC2/DC2/npv100/rcpspdc552.rcp",
    "./data/DC2/DC2/npv100/rcpspdc588.rcp",
    "./data/DC2/DC2/npv100/rcpspdc670.rcp",
    "./data/DC2/DC2/npv100/rcpspdc706.rcp",
]

# Dla każdej z tych instancji pokazać jak rozwiązano je w literaturze
# I jak my rozwiązywaliśmy, zestawić wyniki, rozwiązania, niepoprawne itd.
#
# do poniedziałku rano:
# jakimi metodami i kiedy rozwiązano te instancje z DC2
# ile shotów, optymalnych i tak dalej

# na koniec
# zapisać jakie algorytmy są użwywane do dzielenia problemów na podproblemy w HSS
# sprawdzić czy jest opisane jak ustawiają chain strength i tak dalej

limits = [30, 70, 90]


# if __name__ == '__main__':
#     args = parse_args()
#
#     with open(args.stats_file, 'w') as file:
#         writer = csv.writer(file)
#         writer.writerow("instance,num_of_tasks,num_of_resources,max_length,result,limit")
#         for inst in tqdm(instances):
#             for limit in limits:
#                 instance = read_rcp_instance(inst)
#                 instance.calculate_bounds()
#                 # print(instance)
#                 cqm = instance.get_cqm()
#                 sampler = LeapHybridCQMSampler()
#                 sampleset = sampler.sample_cqm(cqm, time_limit=limit, label=f'PSP name={instance.name}')
#                 data = (
#                     instance.name,
#                     str(len(instance.tasks)),
#                     str(len(instance.resource_caps)),
#                     str(instance.max_length),
#                     str(sampleset.first.sample[f"t{instance.tasks[-1].number}"]),
#                     str(limit),
#                 )
#                 writer.writerow(",".join(data))
#                 # print('Final makespan:', sampleset.first.sample[f"t{instance.tasks[-1].number}"])
#                 # print('Is feasible:', cqm.check_feasible(sampleset.first.sample))

if __name__ == "__main__":
    args = parse_args()
    instance = read_rcp_instance(args.instance)
    instance.calculate_bounds()
    from pprint import pprint

    print(len(instance.__dict__))
    cqm = instance.get_cqm()
    sampler = LeapHybridCQMSampler()
    sampleset = sampler.sample_cqm(
        cqm,
        time_limit=args.time_limit,
        label=f"PSP name={instance.name}, limit={args.time_limit}",
    )
    print(sampleset.first.sample)
    print("Final makespan:", sampleset.first.sample[f"t{instance.tasks[-1].number}"])
    print("Is feasible:", cqm.check_feasible(sampleset.first.sample))

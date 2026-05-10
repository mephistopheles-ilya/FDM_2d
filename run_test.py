import subprocess
import functools
import time

print = functools.partial(print, flush=True)


def print_table_nested (pp, mu, table, h, func, NP):
        print("\\begin{tabular}{", end="")
        for _ in range(2):
            print("|l", end="")
        print("|}")
        print("\\hline")
        if pp == 1.4:
            print("\\multicolumn{", 2 ,"}{|c|}{$\\mu =", str(mu), " \\quad  p(\\rho) = \\rho^{1.4} \\quad ", func, " \\quad \\tau = h = ", str (h),  "$}\\\\", sep="")
        else:
            print("\\multicolumn{", 2 ,"}{|c|}{$\\mu =", str(mu), " \\quad  p(\\rho) =", str(pp), " \\rho \\quad ", func, " \\quad\\tau = h = ", str (h), "$}\\\\", sep="")
        print("\\hline")

        for k in range (1, NP + 1):
            print("$ np = ", str(k) ,"$ ", sep="", end="")
            for norm in range(3):
                print(" & $", str(table.get((h, func, k), ["nan", "nan", "nan"])[norm]), "$\\\\", sep="")
            print("\\hline")
        print("$ real diff ", "$ ", sep="", end="")
        for norm in range(3):
            print(" & $", str(table.get((h, func, 0), ["nan", "nan", "nan"])[norm]), "$\\\\", sep="")
        print("\\hline")
        print("\\end{tabular}")
        print()

def print_table (pp, mu, table, H, HT, func):
        print("\\begin{tabular}{", end="")
        for _ in range(len(H)+ 1):
            print("|l", end="")
        print("|}")
        print("\\hline")
        if pp == 1.4:
            print("\\multicolumn{", len(H) + 1,"}{|c|}{$\\mu =", str(mu), " \\quad  p(\\rho) = \\rho^{1.4} \\quad ", func, "$}\\\\", sep="")
        else:
            print("\\multicolumn{", len(H) + 1 ,"}{|c|}{$\\mu =", str(mu), " \\quad  p(\\rho) =", str(pp), " \\rho \\quad ", func, "$}\\\\", sep="")
        print("\\hline")
        print("$\\tau\\setminus h$", end="")
        for h in H:
            print(" & $", str(h) ,"$", sep="", end="")
        print("\\\\")
        print("\\hline")

        for ht in HT:
            print("$", str(ht) ,"$ ", sep="", end="")
            for norm in range(3):
                for h in H:
                    print(" & $", str(table.get((h, ht, func), ["nan", "nan", "nan"])[norm]), "$", end="", sep="")
                print("\\\\")
            print("\\hline")
        print("\\end{tabular}")
        print()

HT = [0.025, 0.0125, 0.00625, 0.003125, 0.0015625]
H = [0.025, 0.0125, 0.00625, 0.003125, 0.0015625]
MU = [0.1, 0.01, 0.001]
PP = [1, 10, 100, 1.4]

solver = 0


start_time = time.time ()
for pp in PP:
    for mu in MU:
        table = {}
        for h in H:
            for ht in HT:
                result = subprocess.run([
                    "./a.out", 
                    "--hx=" + str(h), 
                    "--hy=" + str(h), 
                    "--ht=" + str(ht),
                    "--mu=" + str(mu), 
                    "--pp=" + str(pp),
                    "--eps=1e-8",
                    "--maxit=2000",
                    "--solver=" + str(solver),
                    "--np=0"
                ], capture_output=True, text=True)

                for line in result.stdout.split('\n'):
                    words = line.split()
                    if "C_norm" in words:
                        table[(h, ht, "G")] = [words[3]]
                        table[(h, ht, "V1")] = [words[6]]
                        table[(h, ht, "V2")] = [words[9]]
                    if "L2_norm" in words:
                        table[(h, ht, "G")].append (words[3])
                        table[(h, ht, "V1")].append (words[6])
                        table[(h, ht, "V2")].append (words[9])
                    if "W1_norm" in words:
                        table[(h, ht, "G")].append (words[3])
                        table[(h, ht, "V1")].append (words[6])
                        table[(h, ht, "V2")].append (words[9])
        print_table (pp, mu, table, H, HT, "G")
        print_table (pp, mu, table, H, HT, "V1")
        print_table (pp, mu, table, H, HT, "V2")

H = [0.0125]
MU = [0.1]
PP = [1, 10, 1.4]

solver = 0

NP = 3

for pp in PP:
    for mu in MU:
        table = {}
        for h in H:
            result = subprocess.run([
                "./a.out", 
                "--hx=" + str(h), 
                "--hy=" + str(h), 
                "--ht=" + str(h),
                "--mu=" + str(mu), 
                "--pp=" + str(pp),
                "--eps=1e-8",
                "--maxit=2000",
                "--solver=" + str(solver),
                "--np=" + str(NP)
            ], capture_output=True, text=True)

            for line in result.stdout.split('\n'):
                words = line.split()
                if "C_norm" in words:
                    table[(h, "G", 0)] = [words[3]]
                    table[(h, "V1", 0)] = [words[6]]
                    table[(h, "V2", 0)] = [words[9]]
                if "L2_norm" in words:
                    table[(h, "G", 0)].append (words[3])
                    table[(h, "V1", 0)].append (words[6])
                    table[(h, "V2", 0)].append (words[9])
                if "W1_norm" in words:
                    table[(h, "G", 0)].append (words[3])
                    table[(h, "V1", 0)].append (words[6])
                    table[(h, "V2", 0)].append (words[9])
                for k in range (1, NP + 1):
                    if "C_norm" + str(k) in words:
                        table[(h, "G", k)] = [words[3]]
                        table[(h, "V1", k)] = [words[6]]
                        table[(h, "V2", k)] = [words[9]]
                    if "L2_norm" + str(k) in words:
                        table[(h, "G", k)].append (words[3])
                        table[(h, "V1", k)].append (words[6])
                        table[(h, "V2", k)].append (words[9])
                    if "W1_norm" + str(k) in words:
                        table[(h, "G", k)].append (words[3])
                        table[(h, "V1", k)].append (words[6])
                        table[(h, "V2", k)].append (words[9])
            print_table_nested (pp, mu, table, h, "G", NP)
            print_table_nested (pp, mu, table, h, "V1", NP)
            print_table_nested (pp, mu, table, h, "V2", NP)

end_time = time.time ()
execution_time = end_time - start_time

print(f"Elapsed: {execution_time} seconds")



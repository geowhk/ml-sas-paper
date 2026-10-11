# Runtime and resource guidance

Times below distinguish original-run records from repository-validation measurements from September–October 2026.
They are not guaranteed completion times. CPU count, memory, operating system, package
builds and contention affect elapsed time; increasing workers does not imply linear speedup.

| Work | Measurement / expectation |
|---|---|
| Saved end-to-end reproduction | 22.8 seconds measured, dependencies already installed; includes F1 bootstrap, Appendix C and ten figures |
| Repeated DGP 2–4, 100 datasets each | Historical batch: 2026-07-20 06:28:07 UTC to 2026-07-21 02:00:28 UTC = 19 h 32 min 21 s; historical parallelism/hardware not fully recorded here |
| Empirical base | Historical run: 05:14:28 to 05:15:53 UTC on 2026-07-15 = 85 s, configuration requested 48 workers |
| Empirical conditional inference | Historical run: 03:54:03 to 03:54:11 UTC on 2026-07-22 = 8 s |
| First dependency install | Separate cost; downloads and source compilation can take many minutes and require system tools |
| Sensitivity / representative DGP | Potentially long; runtime must be measured on the execution machine, not inferred from the saved-result replay |

The current validation host has 10 logical CPUs and 24 GiB RAM. An observed single repeated
DGP process used about 1 GiB RSS early in execution. This is not a guaranteed peak-memory
bound. Empirical 337×99,999 integer permutation indices alone require about 135 MB;
double local-permutation arrays require about 270 MB each, before copies and other objects.

Start with `--workers 1 --parallelism 1`. Increase repeated parallelism only with spare
memory. The full command logs actual stage durations to `reproduction.json` and writes
per-analysis timing rows to `timing.jsonl`; these are the best estimates for a later run
on the same machine. Long stages need not write progress every second.

Full recomputation does not reduce seeds, repetitions, tree counts or permutations to
finish faster. The saved mode is the quick installation/result check. Closing or sleeping
the computer can interrupt a long run: preserve the output folder and use `--resume`.
The repeated stage can resume at replication boundaries. Sensitivity and path/terminal
stages resume at completed-condition boundaries; other incomplete stages restart with
prior incomplete outputs retained. Inspect logs before restarting repeated failures.

## New empirical measurement

Fresh empirical reproduction completed with one worker on the current host:
base analysis 705.2 s, global 99,999-permutation inference 148.4 s, conditional local
99,999-permutation inference 11.4 s, total **865.0 s (14 min 25 s)**. All eight final
CSV tables matched v6 numerically (maximum observed difference zero). Original tree
accumulation groups were retained when reducing executing workers. Other tasks were
running on the host, so these timings include that contention.

The isolated-library saved reproduction took about 31.4 s; the earlier 22.8 s measurement
used the previously installed host packages. Both are measured runs, not performance bounds.

## New repeated-DGP measurement and long sensitivity stages

A fresh DGP 3 replication 91, with both methods and all 99,999 permutations, completed
in **1,651.97 s (27 min 32 s)**. All metric columns matched v6 within 1e-12, including
zero F1 for the path method's no-discovery outcome. macOS `time -l` reported maximum
resident set size 1,715,798,016 bytes (about 1.60 GiB) and peak memory footprint
1,105,594,288 bytes; these are different OS measures, not a memory guarantee.

This one replication cannot predict the whole batch: DGPs and model workloads differ.
The **19 h 32 min historical time covers only the repeated batch**, not the whole full
command. Do not interpret it as an overall completion estimate.

- Hyperparameter sensitivity: 4,500 conditions and 832,500 trees in total.
- Path/terminal sensitivity: 1,800 conditions and 180,000 trees in total.
- Representative DGP and existing-statistics stages add further model and inference work.

The complete sensitivity grids have not been timed again during repository validation. With one
worker, budget potentially multiple days for the full command. No reliable total ETA is
available yet. Each completed sensitivity condition is retained in `conditions/`; an
interrupted condition is recomputed on resume. This can use substantial disk space because
local results are saved per condition. The paper-sized configuration is never reduced by
`--workers` or by resume mode.

The final saved-mode regression after the resume changes completed in 22.2 s
in the isolated R-library environment.

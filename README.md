# TSBProxy
Repository associated to paper "Proxying Betweenness Centrality in Temporal Networks: A comparative Analysis"

# Network file format

The temporal network file must contain one line for each temporal edge (a temporal edge is a triple `(u,v,t)`, where `u` and `v` are two nodes and `t` is one time in which the edge from `u` to `v` is available). All networks are considered as directed: hence, if the graph is undirected, then both the temporal edge `(u,v,t)` and the temporal edge `(v,u,t)` have to be included in the file. Nodes can be identified by any string (not necessarily a number), while times have to be integer values. The file *cannot* include duplicate lines (that is, two identical temporal edges). The three elements of a temporal edge can be separated by any string, which can be specified while using the `load_temporal_graph` function (see below).

# How to use the software

We assume that the `julia` REPL has been started within the directory `www` and the following command has been already executed.

```
include("tsb.jl");
```

## Reading the temporal graph

Assuming that the file `animal.txt` is contained in the directory `graphs` included in the directory `www`, the corresponding temporal network can be loaded by executing the following command (note that, in this case, the separator between the three elements of a temporal edge is the space character).

```
tg = load_temporal_graph("graphs/animal.txt", " ");
```

We can obtain some basic statistics about the temporal network by executing the following command.

```
print_stats(tg, graph_name="Animal interaction network");
```

The result in this case should be the following output.

```
====================================================
Temporal network: Animal interaction network
====================================================
Number of nodes 445
Number temporal of edges 2669
Number of unique time stamps 23
====================================================
```

## Computing the temporal shortest betweenness

The values of the temporal shortest betweenness (in short TSB) of the graph `tg` can be computed by executing the following command.

```
tsb, t = temporal_shortest_betweenness(tg, 1);
```

The second parameter specifies after how many processed nodes a message has to be printed on the console (in order to verify the status of the computation). If this parameter is `0`, then no ouptut is produced. The execution of the above command should require less than one second/. The values returned are the array of the TSB values and the execution time.

The values of the temporal shortest betweenness and teh excution time can be saved in the `scores` and the `times` directory, respectively, (included in the `www` directory) as follows.

```
save_results("animal", "tsb", tsb, t);
```

## Computing the proxies

Similarly to the computation of the TSB values, we can compute and save the values of the four proxies described in the paper, that is, the prefix foremost betweenness (in short, PREFIX), the temporal shortest ego-betweenness (in short, EGOTSB), the prefix foremost ego-betweenness (in short, EGOPREFIX), and the pass-through degree (in short, PTD).

```
prefix, t = temporal_prefix_foremost_betweenness(tg, 100);
save_results("animal", "prefix", prefix, t);
egotsb, t = temporal_ego_betweenness_centrality(tg, 10.0, 100);
save_results("animal", "egotsb", egotsb, t);
egoprefix, t = temporal_ego_prefix_foremost_betweenness(tg, 10.0, 100);
save_results("animal", "egoprefix", egoprefix, t);
ptd, t = pass_through_degree(tg);
save_results("animal", "ptd", ptd, t);
```

## Computing the ONBRA approximation

In the case of ONBRA, we have to decide the size of the sample of node pairs to be used to apply the sampling approximation method. In the paper, we consider three possible sample sizes which should made the execution time of the ONBRA values equal to (respectively, half of and twice) the executon time of PREFIX. To this aim, we first compute the average execution time of ONBRA with a sample of one pair.

```
onbra_v, t = onbra(tg, 10, 0);
```

The average execution time of ONBRA with a sample of one pair is given by `t[1]`. We then divide the execution time of PREFIX by this value in order to obtain the sample size which would cause an execution time approximately equal to PREFIX.

```
ss = Int64(round(read_time("animal", "prefix", -1.0)/t[1];digits=0));
```

We can now compute the ONBRA values with a sample size equal to `ss`, ``2*ss``, and ``0.5*ss``, respectively, and save the values and the times (for each sample_size, we perform 10 experiments).

```
for e in 1:10 onbra_v, t = onbra(tg, ss, 0); save_onbra_results("animal", onbra_v, "equal", ss, e, t); end
ss = 2*ss; for e in 1:10 onbra_v, t = onbra(tg, ss, 0); save_onbra_results("animal", onbra_v, "twice", ss, e, t); end
ss = div(ss, 4); for e in 1:10 onbra_v, t = onbra(tg, ss, 0); save_onbra_results("animal", onbra_v, "half", ss, e, t); end
```

## Analysing the quality of the proxies and of ONBRA

We are now ready to compare the four proxies and ONBRA both in terms of the execution times and in terms of their quality as a proxy. To this last aim, we use both the Spearman correlation index between the different rankings and the minimum h such that the first k nodes in the TSB rankings are includedin the first h position of the ranking of a proxy (we consider h=10, 50, and 100, and h equal to 0.01% of the number of nodes).

```
analyse_all_but_onbra(["animal"]);
analyse_all_onbras(["animal"]);
merge_analysis(["animal"]);
```

After the execution of the above commands, a file `results.csv` is produced in the folder `evaluation/animal`. This file (which is a CSV file with separator `:`) contains all the comparative analysis results. 

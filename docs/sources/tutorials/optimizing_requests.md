# How to Optimize Resource Requests (WIP)

This page outlines some considerations for balancing job performance and
queue wait times.

It's a good idea to become familiar with
[how the clusters are partitioned][nodes-accessible].
This can be useful in maximizing the number of nodes that can potentially
pick up your job. More detailed information about each cluster can be found
on their respective software page (e.g., [Cedar][node-characteristics]).

For monitoring, the live status of nodes on a cluster, there is the
`partition-stats` command:

```bash
$ partition-stats
Node type |                     Max walltime
       |   3 hr   |  12 hr  |  24 hr  |  72 hr  |  168 hr |  672 hr |
----------|-------------------------------------------------------------
       Number of Queued Jobs by partition Type (by node:by core)
----------|-------------------------------------------------------------
Regular   |   12:170 |  69:7066|  70:7335| 386:961 |  59:509 |   5:165 |
Large Mem |    0:0   |   0:0   |   0:0   |   0:15  |   0:1   |   0:4   |
GPU       |    5:14  |   3:8   |  21:1   | 177:110 |   1:5   |   1:1   |
----------|-------------------------------------------------------------
Number of Running Jobs by partition Type (by node:by core)
----------|-------------------------------------------------------------
Regular   |    8:32  |  10:854 |  84:10  |  15:65  |   0:674 |   1:26  |
Large Mem |    0:0   |   0:0   |   0:0   |   0:1   |   0:0   |   0:0   |
GPU       |    5:0   |   2:13  |  47:20  |  19:18  |   0:3   |   0:0   |
----------|-------------------------------------------------------------
       Number of Idle nodes by partition Type (by node:by core)
----------|-------------------------------------------------------------
Regular   |   16:9   |  15:8   |  15:8   |   7:0   |   2:0   |   0:0   |
Large Mem |    3:1   |   3:1   |   0:0   |   0:0   |   0:0   |   0:0   |
GPU       |    0:0   |   0:0   |   0:0   |   0:0   |   0:0   |   0:0   |
----------|-------------------------------------------------------------
       Total Number of nodes by partition Type (by node:by core)
----------|-------------------------------------------------------------
Regular   |  871:431 | 851:411 | 821:391 | 636:276 | 281:164 |  90:50  |
Large Mem |   27:12  |  27:12  |  24:11  |  20:3   |   4:3   |   3:2   |
GPU       |  156:78  | 156:78  | 144:72  | 104:52  |  13:12  |  13:12  |
----------|-------------------------------------------------------------
```

[nodes-accessible]: https://docs.alliancecan.ca/wiki/Job_scheduling_policies#Percentage_of_the_nodes_you_have_access_to
[node-characteristics]: https://docs.alliancecan.ca/wiki/Cedar#Node_characteristics

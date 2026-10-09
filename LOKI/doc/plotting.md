# Plotting

:::{note}
This page is from the earlier LOKI simulation guide and has not been checked since.
:::

The histograms of an analysis (`bcsloki_sans.shist`) can be browsed with

```sh
sb_simplehists_browse bcsloki_sans.shist
```

`sb_loki_plotbcsresults` ({sbpkg}`LOKI/scripts/plotBcsResults`) makes the figures of the LoKI detector reports from
`.shist` files (it needs LaTeX and scipy: `python3 -m pip install scipy`). Without `-pdf` or `-png` it shows the
figure; the numbers of the tabular figures are printed in the terminal. All options: `sb_loki_plotbcsresults -h`.

| Figure | Command |
| --- | --- |
| Polar angle of the neutrons at the source, at conversion and at detection | `sb_loki_plotbcsresults -fn ang3_neutron1e8_isotheta.shist --overlayTheta` |
| Polar angle at detection, per bank | `sb_loki_plotbcsresults -fn ang3_neutron1e8_isotheta.shist` |
| Detection efficiency for different numbers of layers | `sb_loki_plotbcsresults -fn ang1p8_neutron1e8.shist ang3_neutron1e8.shist ang5_neutron1e8.shist ang11_neutron1e8.shist --efficiency` |
| Detection efficiency per bank | `... --bankEfficiency` (same files) |
| Where the neutrons end (absorption), all neutrons | `... --absorption` |
| Where the neutrons end, of those entering the aluminium of a tube | `... --absorption --noAvoid` |

The `.shist` files of the reports are in the
[ESS OwnCloud](https://project.esss.dk/owncloud/index.php/s/DABtPMR4LRWGgIp).

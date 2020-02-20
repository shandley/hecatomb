# Hecatomb databases

This is where the databases probably reside

# Installing the databases

We use several databases to remove contaminants (as shown in the config file section). We provide a compressed tarball that you can easily install yourself, or if you just run snakemake without installing the databases, it will install them for you! (And if they are installed, it will skip that step).

However, we know that sometimes downloading large files over the net causes issues, and so here are the steps to manually install the databases:

1. If you do not want to use the default location, create the directory that you want to use and change into that directory. This is the directory used in `sample_config.yaml` but you probably are not user `redwards` and so change that part as appropriate:

```
cd databases
```

2. Download the databases file. If you have `curl`, I would recommend that:

```
curl -LO https://edwards.sdsu.edu/CERVAID/databases/hecatomb.databases.tar.bz2
```

or if not, try wget

```
wget https://edwards.sdsu.edu/CERVAID/databases/hecatomb.databases.tar.bz2
```

3. Extact the databases:

```
tar xf hecatomb.databases.tar.bz2
```

You're done! That was easy!






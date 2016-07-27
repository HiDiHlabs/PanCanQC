ICGC PanCanQC Docker
===================

This is the first version of the Docker to calculate the quality control values used in the ICGC PanCan project.

Installation
----------------

First of all you need to have a running Docker installation.

Then clone this repository.

```
$ cd PanCanQC/
$ docker build -t "pcawg-qc" .
```

Running the Docker
----------------------------

Then just run the Docker:
```
$ docker run \
		-v /path/to/your/tumor.bam:/home/pcawg/data/bams/tumor.bam \
		-v /path/to/your/control.bam:/home/pcawg/data/bams/control.bam \
		-v /path/to/results_dir/:/home/pcawg/data/results/ \
		-t pcawg-qc
```

You should fine a file with the results in the `/path/to/results_dir/` directory.


Run Segmentation on Dataset with Custom Cellpose Model:
```
vpt --verbose --log-level 1 run-segmentation --segmentation-algorithm cellpose2_custom.json --input-images="images/mosaic_(?P<stain>[\w-]+)_z(?P<z>[0-9]+).tif" --input-micron-to-mosaic images/micron_to_mosaic_pixel_transform.csv --output-path analysis_outputs --tile-size 2400 --tile-overlap 200
```

Create New Cell by Gene Matrix:
```
vpt --verbose partition-transcripts --input-boundaries analysis_outputs/cellpose2_micron_space.parquet --input-transcripts detected_transcripts.csv --output-entity-by-gene analysis_outputs/cell_by_gene.csv
```

Derive Metadata for Resegmented Cells:
```
vpt --verbose derive-entity-metadata --input-boundaries analysis_outputs/cellpose2_micron_space.parquet --output-metadata analysis_outputs/cell_metadata.csv
```

Update .vzg File for the MERSCOPE Visualizer:
```
vpt --verbose --processes 8 update-vzg --input-vzg input_default.vzg --input-boundaries analysis_outputs/cellpose2_micron_space.parquet --input-entity-by-gene analysis_outputs/cell_by_gene.csv --output-vzg analysis_outputs/output_custom.vzg --temp-path temp --overwrite
```

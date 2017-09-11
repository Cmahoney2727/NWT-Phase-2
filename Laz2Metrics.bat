::
:: a batch script for converting raw flight lines (not tiled) into
:: a number of products with a tile-based multi-core batch pipeline
::

::
:: specify parameters
::

:: allows you to run the script from other folders
set PATH=C:\LAStools\bin;

:: here we specify the directory (e.g. the folder) in which the
:: original raw flight lines files are stored
set RAW_FLIGHT_LINES=E:\CFS\ALS_Deliverables\0731a_T03_YZF_YFS

:: here we specify in which format the original raw flight lines 
:: files are stored in the RAW_FLIGHT_LINES folder
set RAW_FORMAT=laz

:: here we specify the directory (e.g. the folder) in which the
:: temporary files are stored
set TEMP_FILES=X:\2016Titan\lastools_temp

:: here we specify the directory (e.g. the folder) in which the
:: resulting output files are stored
set OUTPUT_FILES=E:\CFS\ALS_Deliverables\0731a_T03_YZF_YFS\Outputs

:: here we specify the number of cores we want to run on
set NUM_CORES=4

:: here we specify the target tile size of the tiling
set TILE_SIZE=2000

:: here we specify the target buffer size for each tile
set TILE_BUFFER=1

::rmdir %TEMP_FILES% /s /q

::
:: start processing
::

::goto Counts

:: create buffered tiling from orthometric strips

mkdir %OUTPUT_FILES%\tiles_raw
lastile -i %RAW_FLIGHT_LINES%\*.laz ^
       -tile_size %TILE_SIZE% -buffer %TILE_BUFFER% ^
       -o %OUTPUT_FILES%\tiles_raw\tile -olaz ^
	   -cores %NUM_CORES%

:: ground classify all orthometric tiles

mkdir %OUTPUT_FILES%\tiles_ground_wilderness
lasground -i %OUTPUT_FILES%\tiles_raw\tile*.laz ^
          -wilderness -step 10 ^
          -odir %OUTPUT_FILES%\tiles_ground_wilderness -olaz ^
          -cores %NUM_CORES%

:: remove low and high outliers

mkdir %OUTPUT_FILES%\height_denoised
lasheight -i %OUTPUT_FILES%\tiles_ground_wilderness\*.laz ^
          -replace_z -drop_above 40 -drop_below -0.1 ^
          -odir %OUTPUT_FILES%\height_denoised -olaz ^
          -cores %NUM_CORES%
		  
::Counts

:: building & veggy classify all ground and denoised tiles

mkdir %OUTPUT_FILES%\classified_buffered_tiles
lasclassify -i %OUTPUT_FILES%\height_denoised\*.laz ^
           -odir %OUTPUT_FILES%\classified_buffered_tiles -olaz ^
           -cores %NUM_CORES%

:: Forest Metrics

mkdir %OUTPUT_FILES%\CanopyOutputs
lascanopy -i %OUTPUT_FILES%\classified_buffered_tiles\*.laz ^
			-step 25 -height_cutoff 2 ^
			-avg -std -ske -kur -qav -min -max ^
			-p 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 ^
			-int_avg -int_std -int_ske -int_kur -int_qav -int_min -int_max ^
			-int_p 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 ^
			-use_tile_bb ^
			-odir %OUTPUT_FILES%\CanopyOutputs -otif ^
            -cores %NUM_CORES%
		
:: Count of all returns > 2 m

lascanopy -i %OUTPUT_FILES%\classified_buffered_tiles\*.laz ^
			-step 25 -c 2 40 -use_tile_bb ^
			-odir %OUTPUT_FILES%\CanopyOutputs -otif ^
            -odix _All2m -cores %NUM_CORES%
			
:: Count of first returns > 2 m

lascanopy -i %OUTPUT_FILES%\classified_buffered_tiles\*.laz ^
			-step 25 -first_only -c 2 40 -use_tile_bb ^
			-odir %OUTPUT_FILES%\CanopyOutputs -otif ^
            -odix _First2m -cores %NUM_CORES%

:: Count of all returns

lascanopy -i %OUTPUT_FILES%\classified_buffered_tiles\*.laz ^
			-step 25 -height_cutoff 0 -c 0 40 -use_tile_bb ^
			-odir %OUTPUT_FILES%\CanopyOutputs -otif ^
            -odix _TotalAll -cores %NUM_CORES%

:: Count of first returns

lascanopy -i %OUTPUT_FILES%\classified_buffered_tiles\*.laz ^
			-step 25 -height_cutoff 0 -first_only -c 0 40 -use_tile_bb ^
			-odir %OUTPUT_FILES%\CanopyOutputs -otif ^
            -odix _TotalFirst -cores %NUM_CORES%

::Elevation metrics 

lascanopy -i %OUTPUT_FILES%\tiles_ground_wilderness\*.laz ^
			-step 25 -keep_class 2 -use_tile_bb ^
			-min -max -avg -std ^
			-odir %OUTPUT_FILES%\CanopyOutputs -otif ^
            -odix _Elevation -cores %NUM_CORES%
		  
echo "bye bye"

goto the_end

:the_end

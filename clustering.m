%export LD_PRELOAD=/lib64/libfreetype.so
%matlab -nodisplay -nosplash -nodesktop -r "run('./clustering.m'); exit;" | tail -n +11

% Clustering de espigas obtenidas por los cuadernos de julia
wave_clus_dir = '/home/isabel/Dropbox/Codigo_lab19/m/externals/';

spikes_dir = '/home/isabel/Dropbox/01-06-2022/GDEspigas';
cluster_dir = regexprep( spikes_dir, 'Espigas', 'Times' );

if isfolder( cluster_dir ) == 1;
else
    mkdir( cluster_dir );
end

addpath( genpath( wave_clus_dir ) ); addpath( spikes_dir ); addpath( cluster_dir ); cd( cluster_dir );
chs = dir( fullfile( spikes_dir, '*.mat' ) ); chs = { chs.name }'; chs = strcat( spikes_dir, '/', chs );
errores = { }; j    = 1; i = 1;

while j < length( chs ) && i < length( chs )
    for i = j:length( chs )
        canal = chs{ i };
        msn = strcat( 'Progreso: ', num2str( ceil( ( i/length( chs ) )*100 ) ), ' %' );
        try
            Do_clustering( canal ); disp( msn );
        catch e %#ok<NASGU>
            errores = cat( 1, errores, i ); j = i + 1;
        end
    end
end

save( 'erroresGD.mat', 'errores' );


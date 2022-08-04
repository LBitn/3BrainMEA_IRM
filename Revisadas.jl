"""
    Coleccion de funciones para analisis de archivos BRW y BXR generados con el programa BrainWave de la compañia 3Brain.
	Laboratorio 19 del CINVESTAV a cargo del Dr. Rafael Gutierrez Aguilar.
	Trabajo desarrollado principalmente por Isabel Romero-Maldonado (2020 - )
	isabelrm.biofisica@gmail.com
	https://github.com/LBitn
"""

module Revisadas

# ] add HDF5 JLD DelimitedFiles Statistics Plots HistogramThresholding StatsBase ImageContrastAdjustment Suppressor ImageIO Distributions DSP MAT Dates LinearAlgebra DataFrames

using HDF5
using JLD
using DelimitedFiles
using Statistics
using StatsBase
using Distributions
using Plots
using HistogramThresholding
using ImageContrastAdjustment
using Suppressor
using ImageIO
using DSP
using MAT
using Dates
using LinearAlgebra

# ------------------------------------------------------------------------------------------------------------ #

export minus, searchdir, ms2frames, frames2sec, frames2msec, indexFrCh, matchrow

export div_ab, tup2arr, A_minus_B, IntfromSting, Index2Coordinates, Z0

export BRWDescription, Digital2Analogue

export cortar_cachos, tamaño_cachos, neighborgs
export reparar_rotacion, remover_canales, reparar_saturacion

export ΔV, Thresholding, Figures4
export obtener_grupos, reverberacion

export cortar_Matriz, FiguraGrupos, Figuras_FFinal9

export Parametros, EventosXCanal, EST
export MUA_remez, reparar_separacion
export cortar_espigas, centrar_espigas
export par, MatSpikes, SpikesMat

export meanspike

export Zplot
# ------------------------------------------------------------------------------------------------------------ #
matchrow( a, B ) = findfirst( i -> all( j -> a[ j ] == B[ i, j ], 1:size( B, 2 ) ), 1:size( B, 1 ) );
minus( indx, x ) = setdiff( 1:length( x ), indx );
indexFrCh( total::Int64, Fr::Int64, Ch::Int64 ) = ( total * ( Fr - 1 ) + Ch );
searchdir( path::String, key::String ) = filter( x -> endswith( x, key ), readdir( path; join = true ) );
frames2sec( frames::Real, SamplingRate::Real ) = floor( frames / SamplingRate );
frames2msec( frames::Real, SamplingRate::Real ) = ceil( ( frames / SamplingRate ) * 1000 );
# ------------------------------------------------------------------------------------------------------------ #

"""
"""
function ms2frames( time::Real, SamplingRate::Real )
    if time  != 0
        x = ceil( Int, ( time * SamplingRate ) / 1000 );
    else
        x = 1
    end
    return x
end

"""
    BRWDescription( FILEBRW::String ) -> Dict{ String, Any }
        Lee la version de BRW y si es la última ( 102 ), extrae un diccionario con las variables de interés
        También crea una carpeta Info donde guarda un archivo .jld con ese diccionario para futuros usos.
"""
function BRWDescription( FILEBRW::String )

    BRW = h5open( FILEBRW, "r" );
    Description = read( open_attribute( BRW, "Description" ) );
    BRWversion = read_attribute( open_group( BRW, "3BData" ), "Version" );
    BRWsize = ( ( stat( FILEBRW ).size ) / 1000000 ) / 1024;

    if BRWversion == 102
        RecVars = read( open_group(  BRW, "3BRecInfo" ), "3BRecVars" );
        MeaChip = read( open_group(  BRW, "3BRecInfo" ), "3BMeaChip" );
        MeaChip[ "NRows" ] = convert( Int, MeaChip[ "NRows" ][ 1 ] );
        MeaChip[ "NCols" ] = convert( Int, MeaChip[ "NCols" ][ 1 ] );
        MeaChip[ "Layout" ] =
            reshape( collect( 1:MeaChip[ "NRows" ] * MeaChip[ "NCols" ] ), MeaChip[ "NRows" ], MeaChip[ "NCols" ] );
        Noise = read( open_group( open_group(  BRW, "3BData" ), "3BInfo"), "3BNoise" );
        RawPath = Dict( "Raw" => "3BData/Raw" );
        filename = h5open( FILEBRW, "r" ).filename;
        Description = Dict( "Description" => Description, "Version" => BRWversion, "BRW" => filename );
        SignalInversion = RecVars[ "SignalInversion" ][ ];
        MaxVolt = RecVars[ "MaxVolt" ][ ]; MinVolt = RecVars[ "MinVolt" ][ ]; BitDepth = RecVars[ "BitDepth" ][ ];
        ADCCountsToMV = SignalInversion * ( ( MaxVolt - MinVolt )/ 2^BitDepth );
        MVOffset = SignalInversion * MinVolt;
        Extras = Dict( "MVOffset" => MVOffset, "ADCCountsToMV" => ADCCountsToMV, "BRWsizeGB" => BRWsize );

        # Removing the extra
        delete!( RecVars, "ExperimentType"); delete!( MeaChip, "ROIs" ); delete!( MeaChip, "SysChs" );
        delete!( MeaChip, "MeaType" );

        Variables = merge( Description, RecVars, MeaChip, Noise, RawPath, Extras );

        PATHMAIN = replace( dirname( FILEBRW ), splitpath( FILEBRW )[ end - 1 ] => split( basename( FILEBRW ), "." )[ 1 ] );
        PATHInfo = joinpath( PATHMAIN, "Info" ); mkpath( PATHInfo );
        FILEVARS = joinpath( PATHInfo, "variablesBRW.jld");

        save( FILEVARS, "Variables", Variables );

        close( BRW )
        cd( PATHMAIN )

        println( "You are now working on the new main path: ", PATHMAIN );
        println( "with the file: ", basename( BRW.filename ) );
        println( Description[ "Description" ] );
        println( "HDF5 file size: $BRWsize GB" );

        TEXTVars = copy( Variables );
        delete!( TEXTVars, "StdMean" ); delete!( TEXTVars, "ValidChs" ); delete!( TEXTVars, "Layout" );
        FILEVARS = joinpath( PATHInfo, "variablesBRW.txt");
        writedlm( FILEVARS, TEXTVars );

        return Variables
    else
        println( "You need another function for the: ", BRWversion, " version files, please consult the BRW manual" );
    end
end

"""
"""
function cortar_cachos( Variables::Dict{ String, Any }, Salida_Cachos::String )

    NRecFrames = Variables[ "NRecFrames" ][ 1 ];
    SamplingRate = Variables[ "SamplingRate" ][ 1 ];

    nChs = length( Variables[ "Layout" ] ); Σ = Variables[ "Raw" ];

    Σ = h5open( Variables[ "BRW" ], "r" )[ Σ ];

    ε, ω = tamaño_cachos( NRecFrames, SamplingRate );

    # numero de espacios que ocupa el numero de segundos de duracion del experimento
    n_char = length( string( Int( floor( NRecFrames / SamplingRate ) ) ) );
    # numero de bines
    n = size( ε, 1 );
    # numero de espacios que ocupa el numero de cachos
    n_bins = length( string( n ) );
    #
    for i = 1:n # numero de B a cortar ( 1->4096, 1->ω )
        BIN = Array{ UInt16 }( undef, nChs, ω ); # preallocate
        #= valores correspondientes al BIN especifico.
        El canal 1,1 tiene el frame 1, 4097, 8193...etc =#
        β = collect( ( ε[ i, 1 ] - 1 ):ε[ i, 2 ] );
        for j = 1:ω
            # saca esos frames del machote seguido Σ, y ponlos en array en BIN
            BIN[ :, j ] = Σ[ ( β[ j ] * nChs ) + 1:( nChs * β[ j + 1 ] ) ];
        end
        ini = lpad( string( Int(floor( ( β[ 1 ] )/SamplingRate ) ), "s" ), n_char + 1, "0" );
        eni = lpad( string( Int( floor( ( β[ end ] + 1 ) / SamplingRate ) ), "s" ), n_char + 1, "0" );
        bin_time = string( ini, "-", eni );
        BINname = string( "BIN", lpad( i, n_bins, "0" ), "_", bin_time, ".jld" );
        BINname = joinpath( Salida_Cachos, BINname );
        save( BINname, "data", BIN );
#         println(
#         "Se cortó un segmento de los datos del canal 1 hasta el canal ", nChs,
#         " ( Total: ", nChs, " canales )" );
#         println( " y desde el frame ", ( ε[ i, 1 ] - 1 ), " hasta el frame ", ε[ i, 2 ],
#             "( Total: ", ε[ i, 2 ] - ( ε[ i, 1 ] - 1 ), " frames )" );
#         println( "Con un total de: ", length( BIN ), " puntos" );
    end            
    close( Σ )
end

"""
"""
function tamaño_cachos( NRecFrames, SamplingRate )
    if isinteger( SamplingRate ) # así, se tienen cachos de 1 segundo
        if isinteger( NRecFrames/SamplingRate )
            n = Int.( NRecFrames/SamplingRate );
        end
    elseif isinteger( NRecFrames/floor( SamplingRate ) )
        n = Int.( NRecFrames/floor( SamplingRate ) );
    elseif isinteger( NRecFrames/ceil( SamplingRate ) )
        n = Int.( NRecFrames/ceil( SamplingRate ) );
    else # si no, se hace un merequetengue
        div_T = div_ab( NRecFrames );
        div_sec = div_T/SamplingRate;
        # rango de numero cachos que se quieren cortar, normalmente alto para trabajar a gusto
        hi = 4; lo = 2; # segundos
        if !isempty( div_T )
            # busca uno de los divisores de frames dentro del rango
            selected_divs = div_T[ findall( hi .>= div_sec .>= lo ) ];
        end
        if !isempty( selected_divs ) # si hubo, agarra el primero
            n = Int( NRecFrames/selected_divs[ 1 ] );
        else # si no hubo, uno predeterminado ya que
            n = 60;
        end
    end
    #
    println( " Los ", ( n ), " cachos serán de ", ( ( NRecFrames / SamplingRate ) / n ) ," segundos. " );
    ω = ceil( Int, ( NRecFrames / n ) ); # numero de frames finales (tamaño del cacho en frames)
    ε = Array{ Int64 }( undef, n, 2 ); # preallocate
    ε[ :, 1 ] = collect( 1:ω:NRecFrames ); # inicio y
    ε[ :, 2 ] = ε[ :, 1 ] .+ ω .- 1; # fin en frames de cada cacho (para cortar)
    if !isinteger( NRecFrames / n )
        println( " El ultimo cacho es más chico " )
        ε[ end, 2 ] = NRecFrames; # fin en frames de cada cacho (para cortar)
    end
    return ε, ω
end

"""
    div_ab( n::Int, lo::Int = 1, hi::Int = n )
        Divisores del numero n entre los valores "lo" and "hi", si no estan definidos entonces se toma de 1 a n
    """
function div_ab( n::Int, lo::Int = 1, hi::Int = n )
    ρ = collect( 1:floor( Int, sqrt( n ) ) ) ; # los numeros de 1 en 1 de la raiz cuadrada
    σ1 = findall( n.%ρ .== 0 ); # divisores de la raiz cuadrada ( residuo = 0 )
    σ2 = Int.( ( n ) ./ ( σ1 ) ); # Sacar los pares ( de 100, 2-50, 10-10, etc.)
    σ = sort( unique( vcat( σ1, σ2 ) ) ); # remover duplicados, concatenar, ordenar
    aux1 = @isdefined lo;
    aux2 = @isdefined hi;
    if aux1 && aux2
        rn = σ[ findall( hi .>= σ .>= lo ) ];
        if isempty( rn )
            println(" there is no divisors of $n between $lo and $hi" )
        else
            return rn
        end
    else
        return σ
    end
end

"""
    """
function reparar_rotacion( BIN::Matrix{Float64} )

    if sum( BIN[ 1, : ] ) != 0

        earth = zeros( typeof( BIN[ 1, 1 ] ), size( BIN, 2 ) );
        realearth = matchrow( earth, BIN )[ 1 ];
        BIN = circshift( BIN, ( -1, realearth ) );

        #println( " the refernece channel was at channel ", realearth )
    else
        #println( " The reference channel is in the correct position (channel 1) " )
    end

    return BIN

end

"""
    tup2arr(x::Vector{Tuple})::Vector{NamedTuple{(:Row, :Col), Tuple{Int16, Int16}}}
        Covert the Tuple Vector of Channel coordinates and creates an array with the following content:
            col1 = numero de canal en 4 digitos
            col2 = filas
            col3 = columnas
            numero de canal = ( columna - 1 )*64 + fila
"""
function tup2arr( x::Vector{NamedTuple{(:Row, :Col), Tuple{Int16, Int16}}} )

    aux = zeros(Int, length( x ), 3 );
    for i = 1:length( x )
        aux[ i, 1 ] = ( x[ i ].Col - 1 )*64 + x[ i ].Row;
        aux[ i, 2:3 ] = collect( x[ i ] );
    end
    aux = sortslices( aux, dims = 1 );

    return aux
end

"""
    A_minus_B( A::Array, a::Array )
        Compara cada fila del array chico a con las filas del array grande A
        y si están repetidas las elimina del array grande A
        Devuelve el array grande sin las filas del chico
"""
function A_minus_B( A::Array, a::Array )
    #=
    Compara cada fila del array chico a con las filas del array grande A
    y si están repetidas las elimina del array grande A
    Devuelve el array grande sin las filas del chico
    =#
    n = size( A, 2 );
    m = size( a, 2 );
    if n == m
        for i = 1:size( a , 1 )
            B = reshape( a[ i, : ],( 1, m ) );
            C = ( A .== B );
            C = sum( C, dims = 2 ) .== size( C, 2 );
            D = collect( 1:size( C, 1 ) );
            co = filter!( e -> e .!= 0, vec( Int.( D.*( 1 .- C ) ) ) ) ;
            A = A[ co, : ];
        end
    else
        println("array chico no tiene mismo numero de columnas que array grande")
    end

    return A
end

"""
    remover_canales( data::Matrix{Float64}, thr::Real )
        Calcula aquellos canales cuya desviacion estandar no se comporta como la media.
        Depende del paquete Statistics y Distributions
"""
function remover_canales( data::Matrix{Float64}, thr::Real )

    desviaciones = std(  data, dims = 2 );

    z = fit( Normal, desviaciones );

    limite_derecha = z.μ + thr*z.σ; limite_izquierda = z.μ - thr*z.σ;

    temp = findall( limite_izquierda .< desviaciones .< limite_derecha );

    canales_buenos = zeros( Int, length( temp ) );

    [ canales_buenos[ i, 1 ] = temp[ i ].I[ 1 ] for i in 1:length( temp ) ];

    No_chs = setdiff( 1:size( data, 1 ), canales_buenos );

    return canales_buenos, No_chs

end

"""
    reparar_saturacion( HIthr::Real, LOthr::Real, porcentaje::Real, data::Array )
        elimina las variaciones de voltaje que sobrepasen los umbrales HIthr y LOthr (en micro volts)
        lo hace promediando aquellos canales que estan en la vecindad, siempre y cuando estos no sobrepasen el "porcentaje" de saturacion
    """
function reparar_saturacion( HIthr::Real, LOthr::Real, porcentaje::Real, data::Matrix{Float64} )

    # Toma los canales y frames con saturaciones y repara con la vecindad los que se pueden reparar
    # y los que no los vuelve 0 ( aquellos que estan mas del "porcentaje" de saturacion )
    # depende de las funciones saturacion y neighborgs, y del paquete StatsBase

    porcentaje = porcentaje / 100;
    NOChs = [ ]; NOFrs = [ ]; # preallocation

    # Remover saturaciones positivas y negativas
    # obtener canales saturados para promediacion (Σ).
    ChFrSat = saturacion( data, HIthr, LOthr );

    # Aquellos saturados durante todo el bin son descartados de la lista (ChsGachos)) #
    Σ = Array{ Int64 }( undef, length( countmap( ChFrSat[ :, 1 ] ) ), 2 ); # preallocate
    Σ[ :, 1 ] = Int.( keys( countmap( ChFrSat[ :, 1 ] ) ) ); # que canales
    Σ[ :, 2 ] = Int.( values(countmap( ChFrSat[ :, 1 ] ) ) ); # cuantas veces

    # mas del 50% de frames saturados, se descarta el canal de la lista de promediación
    ChsGachos = Σ[ Σ[ :, 2 ] .>= Int( floor( porcentaje*size( data, 2 ) ) ), 1 ]; # los gachos
    push!( NOChs, ChsGachos ); # lista de todos los canales gachos

    #= obtener Frames saturados para promediacion (Φ).
    Aquellos saturados durante todo el bin son descartados de la lista (FrsGachos)) =#
    Φ = Array{ Int64 }( undef, length( countmap( ChFrSat[ :, 2 ] ) ), 2 ); # preallocate
    Φ[ :, 1 ] = Int.( keys( countmap( ChFrSat[ :, 2 ] ) ) ); # que Frames
    Φ[ :, 2 ] = Int.( values( countmap( ChFrSat[ :, 2 ] ) ) ); # cuantas veces

    # mas del 50% de canales saturados en ese frames (FrsGachos)...son gachos
    FrsGachos = Φ[ Φ[ :, 2 ] .>= Int( floor( porcentaje*size( data, 1 ) ) ), 1 ];
    push!( NOFrs, FrsGachos ); # lista de todos los frames gachos

    #
    # Aquí se quitan los gachos de la lista de reparables
    # final list of channels
    ChFrSat = ChFrSat[ Bool.( 1 .- in.( ChFrSat[ :, 1 ], [ ChsGachos ] ) ), : ];

    # final list of Frames
    ChFrSat = ChFrSat[ Bool.( 1 .- in.( ChFrSat[ :, 2 ], [ FrsGachos ] ) ), : ];

    # ahora ChFrSat contiene solo los canales y frames saturados sin los gachos.
    # Osea, los que se pueden reparar
    for l = 1:size( ChFrSat, 1 )
        Ch = ChFrSat[ l, 1 ]; # channel and
        Fr = ChFrSat[ l, 2 ]; # frame for correction
        # Vecinos del canal gacho
        _, NeighChs = neighborgs( Ch, 1 );

        NeighChsFr = Array{ Int64 }( undef, length( NeighChs ), 2 ); # preallocate
        NeighChsFr[ :, 1 ] = NeighChs; # lista de vecinos
        # cada uno en el frame a promediar
        NeighChsFr[ :, 2 ] = repeat( [ Fr ], length( NeighChs ) );
        # Se remueven de la lista de vecinos los canales gachos
        NeighChsFr = NeighChsFr[
            Bool.( 1 .- in.( NeighChsFr[ :, 1 ], [ ChsGachos ] ) ), : ];
        #= Para evitar reparar el (canal,frame) con sus vecinos igual de saturados se
        remueven los voltajes  superiores a los umbrales establecidos de la lista de
        voltajes vecinos para promediacion =#
        if !isempty( NeighChsFr )
            NeighVoltage = data[ NeighChsFr ][ :, 1 ]; # voltejes de la vecindad
            NeighVoltage = NeighVoltage[
                Bool.( 1 .- ( LOthr .<=  NeighVoltage .>= HIthr ) )
                ]; # voltajes de la vecindad dentro de los umbrales
            if size( NeighVoltage, 1 ) >= 3 # minimo numero de vecinos para promediar
                global data[ Ch, Fr ] = mean( NeighVoltage );
            else
                #= Si no hay suficientes vecinos DENTRO del rango con quienes promediar,
                es mejor matarlo, creo. La otra opción sería promediar con los frames
                inmediatos no saturados del mismo canal...pero no estoy segura. =#
                global data[ Ch, Fr ] = 0;
            end
        else
        # Si no hay suficientes vecinos at all con quienes promediar, es mejor matarlo.
            global data[ Ch, Fr ] = 0;
        end
    end
    #
    # todos los canales gachos y todos los frames gachos se vuelven 0
    data[ ChsGachos, : ] .= 0; data[ :, FrsGachos ] .= 0;
    #
    return data, NOChs, NOFrs
end

"""
    saturacion( data::Array, HIthr::Int, LOthr::Int )
        regresa que canal y en que frame hay saturaciones, osea eventos que se salen por arriba y por abajo de los umbrales establecidos
        HIthr y LOthr en micro volts
    """
function saturacion( data::Matrix{Float64}, HIthr::Real, LOthr::Real )

    saturados = vcat( findall( data .>= HIthr ), findall( data .<= LOthr ) );
    if !isempty( saturados )
        saturados = getindex.( saturados, [1 2] )
    end

    return saturados

end

"""
    neighborgs( center::Int, d::Int ) -> A = Array( (d*2)+1,(d*2)+1), v = vec( 2*((d*2)+1) - 1 );
        a partir del canal (center) se toma la d-vecindad
        A = array donde center es el centro y esta en el orden del chip
        v = mismos canales vecinos que A pero en vector y sin contar el centro
"""
function neighborgs( center::Int64, d::Int64 )

    A = reverse( reshape( collect(1:4096), 64, 64 )', dims = 1 );

    x_c = findall( A .== center )[ ][ 2 ]
    y_c = findall( A .== center )[ ][ 1 ]

    aux = [ ( x_c - d ),( x_c + d ), ( y_c - d ), ( y_c + d ) ]
    aux[ aux .< 1 ] .= 1; aux[ aux .> 64 ] .= 64;

    neigh = A[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    neigh2 = vec( neigh )[ vec( neigh ) .!= center ];

    return neigh, neigh2
end

"""
"""
function Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::Matrix{UInt16} )

    MVOffset = Variables[ "MVOffset" ];
    ADCCountsToMV = Variables[ "ADCCountsToMV" ];

    AnalogValue = @. MVOffset + ( DigitalValue * ADCCountsToMV )

    return AnalogValue

end

"""
"""
function Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::UInt16 )
    MVOffset = Variables[ "MVOffset" ]; ADCCountsToMV = Variables[ "ADCCountsToMV" ];
    AnalogValue = MVOffset + ( DigitalValue * ADCCountsToMV );
    return AnalogValue
end

"""
"""
function ΔV( Variables::Dict{String, Any}, BIN::Matrix{Float64}, ΔT::Int64, descartados::Vector{Int64} )

    SamplingRate = Variables[ "SamplingRate" ][ 1 ];
    ΔT = ms2frames( ΔT, SamplingRate );

    STD = vec( std( ( BIN - circshift( BIN, ( 0, ΔT )  ) ), dims = 2 ) );
    STD2 = copy( STD ); STD2[ descartados ] .= 0;

    return STD, STD2

end

"""
"""
function Thresholding( Variables::Dict{String, Any}, STD::Vector{Float64} )
    nChs = length( Variables[ "Layout" ] );
    W = vec( log.( STD ) );

    edges, conteo = HistogramThresholding.build_histogram( W, length( keys( countmap( W ) ) ) );

    t = zeros( 9 ); Todos = zeros( nChs, length( t ) );

    AbstractImageBinarizationAlgorithm = [
        "UnimodalRosin",
        "MinimumIntermodes",
        "Intermodes",
        "MinimumError",
        "Moments",
        "Otsu",
        "Entropy",
        "Balanced",
        "Yen"
    ];

    @suppress begin

        thr = find_threshold( UnimodalRosin( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 1 ] .= 1;
        t[ 1 ] = thr;
        thr = find_threshold( MinimumIntermodes( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 2 ] .= 1;
        t[ 2 ] = thr;
        thr = find_threshold( Intermodes( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 3 ] .= 1;
        t[ 3 ] = thr;
        thr = find_threshold( MinimumError( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 4 ] .= 1;
        t[ 4 ] = thr;
        thr = find_threshold( Moments( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 5 ] .= 1;
        t[ 5 ] = thr;
        thr = find_threshold( Otsu( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 6 ] .= 1;
        t[ 6 ] = thr;
        thr = find_threshold( Entropy( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 7 ] .= 1;
        t[ 7 ] = thr;
        thr = find_threshold( Balanced( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 8 ] .= 1;
        t[ 8 ] = thr;
        thr = find_threshold( Yen( ), conteo[ 1:end ], edges );
        Todos[ W .>= thr, 9 ] .= 1;
        t[ 9 ] = thr;
    end

    return Todos, t

end

"""
"""
function Figures4( Variables::Dict{ String, Any }, Todos, t, minimos_canales, maximos_canales, STD, FILEFoto, FILESegmento, day, STD2 )

    W = log.( copy( STD ) );

    nChs = length( Variables[ "Layout" ] ); channels = 1:nChs;

    MeaChs2ChIDsMatrix = reverse( reshape( channels, 64, 64 )', dims = 1 );

    AbstractImageBinarizationAlgorithm = [
        "UnimodalRosin",
        "MinimumIntermodes",
        "Intermodes",
        "MinimumError",
        "Moments",
        "Otsu",
        "Entropy",
        "Balanced",
        "Yen"
    ];

    FIGUREstd = heatmap(
        log.( STD2[ MeaChs2ChIDsMatrix ] ),
        aspect_ratio = 1,
        c = :RdBu,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        title = split( basename( FILESegmento ), "." )[ 1 ],
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    aux = sum( Todos, dims = 2 );
    FIGURETodas = heatmap(
        aux[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        title = "Todos los metodos, $( length( unique( aux ) ) - 1 ) grupos",
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    thr = minimum( t[ vec( minimos_canales .>= sum( Todos, dims = 1 ) .<= maximos_canales ) ] );
    Nmetodo = findall( t .== thr )[ 1 ];
    metodo = AbstractImageBinarizationAlgorithm[ Nmetodo ];
    aux = zeros( nChs ); aux[ W .>= thr ] .= 1;
    FIGUREseleccionada = heatmap(
        aux[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        title = metodo,
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    aux = sum( Todos[ :, vec( minimos_canales .>= sum( Todos, dims = 1 ) .<= maximos_canales ) ], dims = 2 );
    FIGURErangos = heatmap(
        aux[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        title = "Dentro del rango, $( length( unique( aux ) ) - 1 ) grupos",
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    if isfile( FILEFoto )

        FIGURERebanada = plot( load( FILEFoto ),
        aspect_ratio = 1,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        title = day,
        titlefont = ( 10, "arial" ) );

        FIGUREFinal = plot( FIGUREstd, FIGURERebanada, FIGURErangos, FIGUREseleccionada, wsize = ( 900, 900 ) );

    else

        FIGUREFinal = plot( FIGUREstd, FIGURETodas, FIGURErangos, FIGUREseleccionada, wsize = ( 900, 900 ) );

    end

    return FIGUREFinal, metodo, Nmetodo
end

"""
    Parametros( SamplingRate::Float64 ) => parametros::Dict{String, Int64}

        lee las variables:
        "Icut", "Lcut", "lF", "fac", "w_post","displace", "stdmin", "HF", "window",
        "stokes", "Vmax", "w_pre"
        Definicion de cada una en el codigo
        Convierte a frames todo lo que este en ms
        Depende de la funcion ms2frames
"""
function Parametros( SamplingRate::Float64 )
    parametros = Dict(
        "window" => ms2frames( window, SamplingRate ),
        "displace" => ms2frames( displace, SamplingRate ),
        "stokes" => ms2frames( stokes, SamplingRate ),
        "w_pre" => ms2frames( w_pre, SamplingRate ),
        "w_post" => ms2frames( w_post, SamplingRate ),
        "Icut" => ms2frames( Icut, SamplingRate ),
        "Lcut" => ms2frames( Lcut, SamplingRate ),
        "stdmin" => stdmin,
        "lF" => lF,
        "HF" => HF,
        "fac" => fac,
        "Vmax" => Vmax,
        "NYQ" => NYQ,
        "order" => order
    );

    return parametros

end

"""
    EST( bin::Array, thr::Real, parametros::Dict{String, Int64} )
        Se obtienen los frames y voltajes que sobrepasan el umbral establecido (thr)
        Si hay eventos supraumbral a n frames de distancia, (stokes) se selecciona aquel que
        tenga menor voltaje.
        Depende del paquete StatsBase

    """
function EST( bin::Array, thr::Real, parametros::Dict{String, Int64} )

    stokes = parametros[ "stokes" ];

    ST = findall( bin .<= thr ); # eventos que pasan el umbral

    if !isempty( ST )

        init = 1;

        while init == 1

            separacion = ST[ 2:end ] .- ST[ 1:( end - 1 ) ]; # distancia entre ellos

            cercanos = findall( separacion .<= stokes ) .+ 1; # cuales estan cerca

            if isempty( cercanos )

                init = 0;

            else

                eliminar = zeros( Int, size( cercanos, 1 ) );

                for i = 1:size( cercanos, 1 )

                    # si el primero es menor que el segundo, quita el segundo

                    if isless( ( bin[ ST[ cercanos[ i ] ] ] ),
                            ( bin[ ST[ cercanos[ i ] - 1 ] ] ) )

                        eliminar[ i ] = cercanos[ i ] - 1;

                    else

                        eliminar[ i ] = cercanos[ i ];

                    end

                end

                if size( eliminar, 1 ) > 1

                    ST[ unique( eliminar ) ] .= 0;

                    filter!( x -> x != 0, ST );

                else

                    ST = ST[ Bool.( ST .!= ST[ eliminar[ 1 ] ] ) ];

                end
            end
        end

        index_parcial = zeros( size( ST, 1 ), 1 );

        index_parcial[ :, 1 ] = ST;

    else

        index_parcial = [ ];

    end

    return index_parcial

end

"""
    EventosXCanal( canal::Array, parametros::Dict{String, Int64} )
        Identifica el numero de eventos en una ventana temporal que se desplaza
        Stokes es la distancia minima entre eventos para ser considerados uno solo
        Depende de la funcion EST, BRW y el paquete Statistics
"""
function EventosXCanal( canal::Array, parametros::Dict{String, Int64} )

    index = [ ];

    window = parametros[ "window" ]; displace = parametros[ "displace" ]; stdmin = parametros[ "stdmin" ];

    i = 1

    Inicio = Int( ( ( ( i - 1 ) * displace ) + 1 ) ); Final = Int( Inicio + window - 1 );
    thrs = 0;

    while Final <= ( length( canal ) - Int( window - 1 ) )

        bin =  canal[ Inicio:Final ];

        if !iszero( bin )

            thr = -1 * stdmin * ( median( abs.( bin ) ) / 0.6745 );

            thrs = vcat( thrs, thr );

            index_parcial = Int.( EST( bin, thr, parametros ) );

            if !isempty( index_parcial )

                index_real = index_parcial .+ Inicio .- 1;

                index = vcat( index, index_real );

            end

        end

        i = i + 1;

        Inicio = Int( ( ( ( i - 1 ) * displace ) + 1 ) ); Final = Int( Inicio + window - 1 );

    end

    thrs = thrs[ 2:end ]; index = unique( index );

    return index, thrs

end

"""
    MUA_remez( canal::Array, parametros::Dict{String, Int64}, SamplingRate::Float64 )
        =>  MUA::Vector{Float64}, canal::Array{Float64, 1}
        toma el canal crudo y lo filtra con un filtro remez con atenuacion de "factor"
        entre las bandas de frecuencia.
        Depende del paquete DSP
    """
function MUA_remez( canal::Array, parametros::Dict{String, Int64}, SamplingRate::Float64 )

    canal = Float64.( canal );

    lF = parametros[ "lF" ]; fac = parametros[ "fac" ]; HF = parametros[ "HF" ];
    order = parametros[ "order" ]; NYQ = parametros[ "NYQ" ];

    bpass = remez(
        ( order + 1 ), [ ( 0, lF - fac ) => 0, ( lF, HF ) => 1, ( HF + fac, NYQ ) => 0 ],
            Hz = SamplingRate );

    MUA = filtfilt( DSP.Filters.PolynomialRatio( bpass, [ 1.0 ] ), canal );

    return MUA, canal
end

"""
    reparar_separacion( index::Vector{Any}, MUA::Vector{Float64}, parametros::Dict{String, Int64} )
        => index_final::Vector{Int64}

        Por algun motivo si la ventana de muestreo es muy corta en el paso previo (EventosXCanal),
        hay redundancia en la deteccion de algunos eventos.
        Este codigo es para volver a aplicar la regla de la distancia minima previamente establecida

"""
function reparar_separacion( index::Vector{Any}, MUA::Vector{Float64}, parametros::Dict{String, Int64} )

    stokes = parametros[ "stokes" ];

    init = 1;

    while init == 1

        separacion = index[ 2:end ] .- index[ 1:( end - 1 ) ]; # distancia entre ellos

        cercanos = findall( separacion .<= stokes ) .+ 1; # cuales estan cerca

        if isempty( cercanos )

            init = 0;

        else

            eliminar = zeros( Int, size( cercanos, 1 ) );

            for i = 1:size( cercanos, 1 )

                # si el primero es menor que el segundo, quita el segundo

                if isless( ( MUA[ index[ cercanos[ i ] ] ] ),
                        ( MUA[ index[ cercanos[ i ] - 1 ] ] ) )

                    eliminar[ i ] = cercanos[ i ] - 1;

                else

                    eliminar[ i ] = cercanos[ i ];

                end

            end

            if size( eliminar, 1 ) > 1

                index[ unique( eliminar ) ] .= 0;

                filter!( x -> x != 0, index );

            else

                index = index[ Bool.( index .!= index[ eliminar[ 1 ] ] ) ];

            end
        end
    end

    index_final = zeros( size( index, 1 ), 1 ); index_final[ :, 1 ] = index;

    return vec( Int.( index_final ) )

end

"""
    cortar_espigas( index_final::Vector{Any}, MUA::Vector{Float64}, parametros::Dict{String, Int64} )
        => spikes::Matrix{Float64}

        toma los centros detectados en el paso EventosXCanal y corta las espigas, con w_pre ms antes
        y w_pos ms despues.

    """
function cortar_espigas( index_final::Vector{Int64}, MUA::Vector{Float64}, parametros::Dict{String, Int64} )

    w_post = parametros[ "w_post" ]; w_pre = parametros[ "w_pre" ];

    index_final = index_final[ Bool.( 1 .- ( ( index_final .+ w_post ) .> length( MUA ) ) ) ];

    index_final = index_final[ Bool.( 1 .- ( ( index_final .- w_pre ) .< 1 ) ) ];

    spikes = zeros( length( index_final ), ( w_pre + w_post + 1) );

    for i = 1:length( index_final )

        spikes[ i, : ] = MUA[ ( index_final .- w_pre )[ i ]:( index_final .+ w_post )[ i ] ];

    end

    return spikes

end

"""
    function centrar_espigas( spikes::Matrix{Float64}, parametros::Dict{String, Int64} )
        => spikes::Matrix{Float64}
        Determina si los centros de las espigas estan movidos unos frames a la izquierda o la derecha y
        los recorre solo para obtener un mejor sorting
        Ignora los extremos
        Depende del paquete StatsBase
    """
function centrar_espigas( spikes::Matrix{Float64}, parametros::Dict{String, Int64} )

    Icut = parametros["Icut"]; Lcut = parametros["Lcut"]; Lcut = size( spikes, 2 ) - Lcut;

    frames_centros = zeros( Int, size( spikes, 1 ) );

    for i = 1:size( spikes, 1 )

        frames_centros[ i ] = findall( minimum( spikes, dims = 2 )[ i ] .== spikes[ i, : ] )[ 1 ];

    end

    moda = StatsBase.mode( frames_centros ); espigas_movidas = findall( frames_centros .!= moda );

    por_corregir = espigas_movidas[ Icut .< frames_centros[ espigas_movidas ] .< Lcut ];

    cuanto_corregir =
        frames_centros[
            espigas_movidas ][ Icut .< frames_centros[ espigas_movidas ] .< Lcut ];

    retrasar = por_corregir[ ( cuanto_corregir .- moda ) .> 0 ];

    adelantar = por_corregir[ ( cuanto_corregir .- moda ) .< 0 ];

    if !isempty( retrasar )

        cuanto_retrasar = frames_centros[ retrasar ] .- moda;

        for i = 1:length( retrasar )

            spikes[
                retrasar[ i ], : ] = circshift(
                    spikes[ retrasar[ i ], : ], ( -1 * cuanto_retrasar[ i ] ) );

        end

    end
    if !isempty( adelantar )

        cuanto_adelantar = frames_centros[ adelantar ] .- moda;

        for i = 1:length( adelantar )

            spikes[
                adelantar[ i ], : ] = circshift(
                    spikes[ adelantar[ i ], : ], ( -1 * cuanto_adelantar[ i ] ) );

        end

    end

    return spikes

end

"""
    par( Variables::Dict{String, Any}, parametros::Dict{String, Int64} )
        => par::Dict{String, Any}
        Obtiene las variables necesarias para rellendar el archivo .mat para hacer el SPC
        Depende del paquete Dates
"""
function par( Variables::Dict{String, Any}, parametros::Dict{String, Int64} )

    right_now =
        Dates.format( round( Dates.DateTime( Dates.now( ) ), Dates.Minute( 15 ) ), Dates.RFC1123Format );

    SamplingRate = Variables[ "SamplingRate" ];

    par = Dict(
            "channels"         => 1,
            "segments_length"  => ( ( ceil( Variables[ "NRecFrames" ][ 1 ] )/ SamplingRate ) ) / 60,
            "sr"               => SamplingRate,
            "tmax"             => "all",
            "tmin"             => 0,
            "w_pre"            => parametros[ "w_pre" ],
            "w_post"           => parametros[ "w_post" ],
            "alignment_window" => parametros[ "w_pre" ] + 1,
            "stdmin"           => parametros[ "stdmin" ],
            "stdmax"           => 35,
            "detect_fmin"      => parametros[ "lF" ],
            "detect_fmax"      => parametros[ "HF" ],
            "sort_fmin"        => parametros[ "lF" ],
            "sort_fmax"        => parametros[ "HF" ],
            "ref_ms"           => 1,
            "detection"        => "neg",
            "int_factor"       => 5,
            "interpolation"    => "y",
            "sort_order"       => 2,
            "detect_order"     => parametros[ "order" ] + 1,
            "detection_date"   => right_now
                );

    return par

end

"""
    MatSpikes( param::Dict{String, Any}, parametros::Dict{Any, Any},
        spikes::Matrix{Float64}, MUA::Vector{Float64} )
        A partir de las variables previamente definidas, crea un archivo .mat para hacer el SPC
        Depende del paquete MAT

"""
function MatSpikes( param::Dict{String, Any}, parametros::Dict{Any, Any},
        spikes::Matrix{Float64}, MUA::Vector{Float64} )

    matwrite( parametros[ "FILEspikes" ],

            Dict(
            "par" => param,
            "threshold" => abs( parametros[ "thr" ] ),
            "index" => reshape(
                    ( ( parametros[ "index" ]./param[ "sr" ][ 1 ] )*1000 ),
                        1, length( parametros[ "index" ] ) ),
            "spikes" => spikes,
            "psegment" => MUA,
            "sr_psegment" => param[ "sr" ][ 1 ],
        ); compress = true

        );

    println( "espigas del canal -> ", parametros[ "Qcanal" ], " listas" );

end

"""
    SpikesMat( Ncanal::Int64, canal::Vector{Float16}, parametros::Dict{String, Int64},
        Variables::Dict{String, Any} )
        Toma un canal crudo, lo filtra, detecta eventos, repara errores en la deteccion, corta los eventos
        los centra, remueve los Outlayers, y crea un archivo .mat para hacer el SPC
        Depende del paquete Statistics
    """
function SpikesMat( Ncanal::Int64, canal::Vector{Float16}, parametros::Dict{String, Int64},
        Variables::Dict{String, Any}, PATHSpikes::String )

    SamplingRate = Variables[ "SamplingRate" ][ 1 ];

    MUA, canal = MUA_remez( canal, parametros, SamplingRate );

    index, thrs = EventosXCanal( MUA, parametros );

    index_final = reparar_separacion( index, MUA, parametros );

    spikes = cortar_espigas( index_final, MUA, parametros ); spikes = centrar_espigas( spikes, parametros );

    # Depuraci'ones
    # Espigas que sobrepasan el limite positivo de voltaje

    Vmax = parametros[ "Vmax" ];

    Outlayers = findall( spikes .> Vmax );

    if !isempty( Outlayers )

        aux = sort( unique( getindex.( Outlayers, [ 1 ] ) ) );

        si_pasan = setdiff( 1:size( spikes, 1 ), aux );

        index_final = index_final[ si_pasan ]; thrs = thrs[ si_pasan ]; spikes = spikes[ si_pasan, : ];

    end

    # Espigas que sobrepasan el limite de desviaci'on estandar

    aux = std( spikes, dims = 2 ); y = ( mean( aux ) + 2*std( aux ) ); si_pasan = findall( vec( aux ) .< y );

    index_final = index_final[ si_pasan ]; thrs = thrs[ si_pasan ]; spikes = spikes[ si_pasan, : ];

    # Valores para el .mat

    thr_final = mean( thrs ); param = par( Variables, parametros );

    FILEspikes = joinpath( PATHSpikes, string( "ch", lpad( Ncanal, 4, "0" ), "_spikes.mat" ) );

    p = Dict(); p[ "FILEspikes" ] = FILEspikes; p[ "index" ] = index_final;

    p[ "thr" ] = thr_final; p[ "Qcanal" ] = Ncanal; parametros = merge( p, parametros );

    MatSpikes( param, parametros, spikes, MUA )

end


"""
"""
function reverberacion( grupos )

    canales_juntos = 0

    for i in grupos

        canales_juntos = vcat( canales_juntos, i )

    end

    canales_juntos = canales_juntos[ canales_juntos .!= 0 ];

    canales_juntos = unique( canales_juntos );

    more_groups = [ ]

    for i in canales_juntos

        ok = 0

        for j = 1:length( grupos )

            if !isempty( findall( grupos[ j ] .== i ) )

                ok = vcat( ok, grupos[ j ] );

            end
        end

        ok = ok[ok .!= 0];

        new_group = sort( unique( ok ) );

        if !isempty( new_group )

            if !isempty( more_groups )

                temp = last( more_groups );

                if !isequal( temp, new_group )

                    push!( more_groups, new_group );
                end

            else

                push!( more_groups, new_group );

            end

        end

    end

    return more_groups

end

"""
"""
function obtener_grupos( seleccion::Vector{Int64} )

    grupos = [ ];

    for i in seleccion

        _, vecinos = neighborgs( i, 1 );

        grupo = sort( intersect( vecinos, seleccion ) );

        if isempty( grupo )

            push!( grupos, i )

        else

            push!( grupos, vcat( grupo, i ) )

        end

    end

    canales_solos = grupos[ findall( length.( grupos ) .== 1 ) ];

    deleteat!( grupos, findall( length.( grupos ) .== 1 ) );

    a = length( grupos ); grupos = reverberacion( grupos ); b = length( grupos );

    while a != b

        a = length( grupos ); grupos = reverberacion( grupos ); b = length( grupos );

    end

    return grupos, canales_solos

end

"""
    IntfromSting( X::String ) -> n
        Busca dentro del texto X si hay numeros que puedan ser convertido a entero.
        Se asume que ser'an un numero continuo de digitos para ser transformados a un entero
        Entonces regresa el entero formado desde el primer hasta el ultimo digito convertible.
"""
function IntfromSting( X::String )

    a = findfirst( isdigit.( collect( X ) ) )
    b = findlast( isdigit.( collect( X ) ) )

    if !isempty( a ) && !isempty( b )
        n = parse( Int, X[ a:b ] )
    else
        println(" There is no integers into this string ")
    end

    return n

end

"""
"""
function Index2Coordinates( x::Int64 )
    nChs = 4096;
    MeaChs2ChIDsMatrix = reverse( reshape( collect( 1:nChs ), 64, 64 )', dims = 1 );

    aux = findall( MeaChs2ChIDsMatrix .== x )[ 1 ]; x1 = aux[ 1 ]; y1 = aux[ 2 ];

    return [ x1, y1 ]
end

"""
"""
function Index2Coordinates( x::Vector{Int64} )
    nChs = 4096;
    MeaChs2ChIDsMatrix = reverse( reshape( collect( 1:nChs ), 64, 64 )', dims = 1 );

    r = zeros( Int, length( x ), 3 )
    for i = 1:length( x )
        aux = findall( MeaChs2ChIDsMatrix .== x[ i ] )[ 1 ]; x1 = aux[ 1 ]; y1 = aux[ 2 ];
        r[ i, : ] = [ x[ i ], x1, y1 ];
    end

    return r
end

# Arreglar estas funciones....
function Figuras_FFinal9( MeaChs2ChIDsMatrix, Todos, FILEFiguraFinal )

    F1 = heatmap(
        Int.( Todos[ :, 1 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );
    F2 = heatmap(
        Int.( Todos[ :, 2 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    F3 = heatmap(
        Int.( Todos[ :, 3 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    F4 = heatmap(
        Int.( Todos[ :, 4 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    F5 = heatmap(
        Int.( Todos[ :, 5 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    F6 = heatmap(
        Int.( Todos[ :, 6 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    F7 = heatmap(
        Int.( Todos[ :, 7 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    F8 = heatmap(
        Int.( Todos[ :, 8 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    F9 = heatmap(
        Int.( Todos[ :, 9 ] )[ MeaChs2ChIDsMatrix ],
        aspect_ratio = 1,
        c = :greys,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        cbar = :none,
        titlefont = ( 10, "arial" ),
    );

    FFinal = plot( F1, F2, F3, F4, F5, F6, F7, F8, F9, layout = ( 3, 3 ), wsize = ( 800, 800 ) );

    savefig( FFinal, FILEFiguraFinal )

end

function FiguraGrupos( FILEFiguraGrupos, Todos, FILEVars )

    Variables = load( FILEVars )[ "Variables" ];

    nChs = length( Variables[ "Layout" ] );

    channels = 1:nChs;

    AbstractImageBinarizationAlgorithm = [
            "UnimodalRosin",
            "MinimumIntermodes",
            "Intermodes",
            "MinimumError",
            "Moments",
            "Otsu",
            "Entropy",
            "Balanced",
            "Yen"
        ];

    NMetodo = findall( maximum( sum( Todos, dims = 1 ) ) .== vec( sum( Todos, dims = 1 ) ) )[ 1 ];

    Metodo = AbstractImageBinarizationAlgorithm[ NMetodo ];

    NCanales_Finales = Int( maximum( sum( Todos, dims = 1 ) ) );

    Canales_Finales = channels[ Bool.( Int.( Todos[ :, NMetodo ] ) ) ];

    grupos, canales_solos = obtener_grupos( Canales_Finales );

    Z = zeros( Int, nChs );

    for i = 1:length( grupos )

        Z[ grupos[ i ] ] .= i;

    end
        Z = reshape( Z, 64, 64 );

        FIGGrupos = heatmap(
        Z',
        aspect_ratio = 1,
        c = :twilight,
        axis = ( [ ], false ),
        wsize = ( 400, 400 ),
        yflip = true,
        cbar = :none

    )

    savefig( FIGGrupos, FILEFiguraGrupos );

    Canales = vcat( grupos...); Canales = vcat( Canales, canales_solos );

    Variables[ "ACD" ] = Canales;
    Variables[ "Metodo" ] = Metodo;
    Variables[ "NCanales_Finales" ] = NCanales_Finales;

    save( FILEVars, "Variables", Variables );

    return Canales, FIGGrupos
end

# -------------------------------------------------------------------------------------------------------- #

function cortar_Matriz( Canales, NRecFrames, FILESNoSat, FILEMatriz )

    Matriz = Array{ Float16 }( undef, length( Canales ), NRecFrames );
    A = 1;
    B = Int( NRecFrames / length( FILESNoSat ) );

    for i = 1:length( FILESNoSat )

        data = load( FILESNoSat[ i ] )[ "data" ];

        Matriz[ :, A:B ] = data[ Canales, : ];

        A = B + 1;
        B = B + Int( NRecFrames/length( FILESNoSat ) );

    end

    # Guardando matriz y que electrodos contiene

    save( FILEMatriz, "Matriz", Matriz );

    mapa_index = hcat( 1:length( Canales ), Canales );

    return Matriz, mapa_index

end

"""
"""
function Z0( X )

    Z = zeros( Int, 4096 );
    Z[ X ] = Z[ X ] .+ 1;
    Z = reverse( reshape( Z, 64, 64 )', dims = 1 );

    return Z
end

"""
"""
function Zplot( Z, cm )

    F = heatmap(
            Z,
            aspect_ratio = 1,
            c = cm,
            axis = ( [ ], false ),
            wsize = ( 400, 400 ),
            cbar = :none );
    return F

end

# -------------------------------------------------------------------------------------------------------- #

function meanspike( spikes::Array, thr::Real )

    x = cor( spikes, dims = 2 ); x[ x .== 1 ] .= 0;

    y = getindex.( findall( UpperTriangular( x ) .>= thr ), [ 1 2 ] );

    CArray = zeros( Float64, size( y, 1 ), size( spikes, 2 ) );

    for i = 1:size( y, 1 )

        CArray[ i, : ] =  mean( [ spikes[ y[ i, : ][ 1 ], : ], spikes[ y[ i, : ][ 2 ], : ] ] );

    end

    mean_spike = vec( mean( CArray, dims = 1 ) );

    return mean_spike

end

# -------------------------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------------------------- #
#
end

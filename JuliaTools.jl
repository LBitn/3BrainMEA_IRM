"""
    Collection of functions for analysis of BRW (HDF5) files generated with the BrainWave program from the company 3Brain.
    Laboratory 19 of the CINVESTAV in charge of Dr. Rafael Gutierrez Aguilar.
    Work developed mainly by Isabel Romero-Maldonado (2020 - )
    isabelrm.biofisica@gmail.com
    https://github.com/LBitn
"""

module JuliaTools

using HDF5
using JLD
using DelimitedFiles
using Statistics
using StatsBase
using Distributions
using Plots

export searchdir
export Sparcity
export Density
export Get_Groups
export neighborgs
export Zplot
export ms2frames
export donoho
export MeanΔxCI
export div_ab
export size_segments
export Get_chunks
export BRWDescription
export Digital2Analogue

"""
    Noise-adaptive Optimal Thresholding
"""
@inline donoho( x ) =  ( median( abs.( x ) ) / 0.6745 );

"""
    searchdir( path::String, key::String ) -> Vector{String}
    Find inside the given path, the files with the key word on their name and returns a vector of strings with the full name of those files (complete path)
"""
@inline searchdir( path::String, key::String ) = filter( x -> endswith( x, key ), readdir( path; join = true ) );

@inline Sparsity( count_non_zeros::Int64, total_elements_of_A::Int64 ) = 1 - count_non_zeros/total_elements_of_A;
@inline Density( W::Vector{Any} ) = length( findall( ( length.( W ) ./ length( W ) ) .>= ( mean( length.( W )./length( W ) ) + 2*std( length.( W )./length( W ) ) ) ) );


"""
    Get_Groups( W::Vector{Int64} ) -> grupos::Vector{Any}, loose::( Vector{Any} )
       Groups of adjacent channels are formed from the initial indexes W. Those channels separated from the bulk of the others is considered loose ones.
"""
function Get_Groups( W::Vector{Int64} )

    grupos = [ ];

    for i in W

        _, vecinos = neighborgs( i, 1 );

        grupo = sort( intersect( vecinos, W ) );

        if isempty( grupo )

            push!( grupos, i )

        else

            push!( grupos, vcat( grupo, i ) )

        end

    end

    loose = grupos[ findall( length.( grupos ) .== 1 ) ];

    deleteat!( grupos, findall( length.( grupos ) .== 1 ) );

    a = length( grupos ); grupos = reverberation( grupos ); b = length( grupos );

    while a != b

        a = length( grupos ); grupos = reverberation( grupos ); b = length( grupos );

    end

    return grupos, loose

end

"""
    reverberation( grupos::Vector{Any} ) -> more_groups::Vector{Any}
        Intermediate step for Get_Groups function
"""
function reverberation( grupos::Vector{Any} )

    adjoining_channels = 0

    for i in grupos

        adjoining_channels = vcat( adjoining_channels, i )

    end

    adjoining_channels = adjoining_channels[ adjoining_channels .!= 0 ];

    adjoining_channels = unique( adjoining_channels );

    more_groups = [ ]

    for i in adjoining_channels

        temporal = 0

        for j = 1:length( grupos )

            if !isempty( findall( grupos[ j ] .== i ) )

                temporal = vcat( temporal, grupos[ j ] );

            end
        end

        temporal = temporal[ temporal .!= 0 ];

        new_group = sort( unique( temporal ) );

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
    neighborgs( C::Int64, d::Int64 ) ->
        -> A = Array( ( d*2 ) + 1, ( d * 2 ) + 1 ), v = vec( 2*( ( d * 2 ) + 1 ) - 1 );
        The d-neighborhood is calculated from the channel (C) as a center
        A = array where C is the center and is in chip order
        v = same neighboring channels as A but in vector form and without C ( 8 channels )
"""
function neighborgs( C::Int64, d::Int64 )

    Layout = reverse( reshape( collect( 1:4096 ), 64, 64 )', dims = 1 );

    x_c = findall( Layout .== C )[ ][ 2 ]; y_c = findall( Layout .== C )[ ][ 1 ];

    aux = [ ( x_c - d ),( x_c + d ), ( y_c - d ), ( y_c + d ) ]
    aux[ aux .< 1 ] .= 1; aux[ aux .> 64 ] .= 64;

    A = Layout[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    v = vec( A )[ vec( A ) .!= C ];

    return A, v
end

"""
    MeanΔxCI( W::Vector, percent::Float64 ) -> xmean, Δx, C1, C2 (::Float64)
    Assumes a Normal distribution. Obtains the confidence interval with "percent" quantile,
    Returns mean, standard error, CI superior and CI inferior
"""
function MeanΔxCI( W::Vector, percent::Float64 )

    xmean, Δx = BinningAnalysis.jackknife( identity, W );
    d = Distributions.Normal( xmean, std( W ) );
    C1 = xmean + ( Δx * Distributions.quantile( d, percent ) );
    C2 = xmean - ( Δx * Distributions.quantile( d, percent ) );

    return xmean, Δx, C1, C2
end

"""
    Zplot( Z::Array, cm::colormap ) -> F::heatmap
        Just a plotting recipe
        cm = colormap (:greys)
        Z = Array(64x64)
"""
# - Auxiliars for ploting
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

"""
    ms2frames( time::Real, SamplingRate::Real ) -> x::Int64
        For conversion from ms to an integer number of frames
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
    Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::UInt16 )
"""
function Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::Matrix{UInt16} )

    MVOffset = Variables[ "MVOffset" ];
    ADCCountsToMV = Variables[ "ADCCountsToMV" ];

    AnalogValue = @. MVOffset + ( DigitalValue * ADCCountsToMV )

    return AnalogValue

end

"""
    Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::UInt16 )
"""
function Digital2Analogue( Variables::Dict{ String, Any }, DigitalValue::UInt16 )
    MVOffset = Variables[ "MVOffset" ]; ADCCountsToMV = Variables[ "ADCCountsToMV" ];
    AnalogValue = MVOffset + ( DigitalValue * ADCCountsToMV );
    return AnalogValue
end

"""
    BRWDescription( FILEBRW::String ) -> Dict{ String, Any }
        Reads the version of BRW and if it is the latest ( 102 ) version, extract a dictionary with the variables of interest
        It also creates an Info folder where it saves a .jld file with that dictionary for future uses.
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
    Get_chunks( Variables::Dict{ String, Any }, Output_Chunks::String ) -> nothing
        Cuts the brw dataset into several more manageable segments
"""
function Get_chunks( Variables::Dict{ String, Any }, Output_Chunks::String )

    NRecFrames = Variables[ "NRecFrames" ][ 1 ];
    SamplingRate = Variables[ "SamplingRate" ][ 1 ];

    nChs = length( Variables[ "Layout" ] ); Σ = Variables[ "Raw" ];

    Σ = h5open( Variables[ "BRW" ], "r" )[ Σ ];

    ε, ω = size_segments( NRecFrames, SamplingRate );

    # number of spaces that occupies the number of seconds of duration of the experiment
    n_char = length( string( Int( floor( NRecFrames / SamplingRate ) ) ) );
    # number of bins
    n = size( ε, 1 );
    # number of spaces occupied by the number of chunks
    n_bins = length( string( n ) );
    #
    for i = 1:n # number of B to cut ( 1->4096, 1->ω )
        BIN = Array{ UInt16 }( undef, nChs, ω ); # preallocate
        #= values corresponding to the specific BIN.
         Channel 1,1 has the frames 1, 4097, 8193...etc =#
        β = collect( ( ε[ i, 1 ] - 1 ):ε[ i, 2 ] );
        for j = 1:ω
            # take those frames out of the vector Σ, and put them in array in BIN
            BIN[ :, j ] = Σ[ ( β[ j ] * nChs ) + 1:( nChs * β[ j + 1 ] ) ];
        end
        ini = lpad( string( Int(floor( ( β[ 1 ] )/SamplingRate ) ), "s" ), n_char + 1, "0" );
        eni = lpad( string( Int( floor( ( β[ end ] + 1 ) / SamplingRate ) ), "s" ), n_char + 1, "0" );
        bin_time = string( ini, "-", eni );
        BINname = string( "BIN", lpad( i, n_bins, "0" ), "_", bin_time, ".jld" );
        BINname = joinpath( Output_Chunks, BINname );
        save( BINname, "data", BIN );
    end
    close( Σ )
end

"""
    size_segments( NRecFrames::Int64, SamplingRate::Float64 ) -> ε::Vector, ω::Inf64
        Determines the ideal size of each chunck

"""
function size_segments( NRecFrames::Int64, SamplingRate::Float64 )
    if isinteger( SamplingRate ) # thus, there are chunks of 1 second
        if isinteger( NRecFrames/SamplingRate )
            n = Int.( NRecFrames/SamplingRate );
        end
    elseif isinteger( NRecFrames/floor( SamplingRate ) )
        n = Int.( NRecFrames/floor( SamplingRate ) );
    elseif isinteger( NRecFrames/ceil( SamplingRate ) )
        n = Int.( NRecFrames/ceil( SamplingRate ) );
    else # if not, a merequetengue is made
        div_T = div_ab( NRecFrames );
        div_sec = div_T/SamplingRate;
        # range of number of chunks to be cut, normally high to work at ease
        hi = 4; lo = 2; # segundos
        if !isempty( div_T )
            # finds one of the n-frames dividers within the range
            selected_divs = div_T[ findall( hi .>= div_sec .>= lo ) ];
        end
        if !isempty( selected_divs ) # if there was, grab the first one
            n = Int( NRecFrames/selected_divs[ 1 ] );
        else # if there wasn't, a default one.
            n = 60;
        end
    end
    #
    println( " The are ", ( n ), " segments, of ", ( ( NRecFrames / SamplingRate ) / n ) ," seconds each." );
    ω = ceil( Int, ( NRecFrames / n ) ); # number of final frames (chunk size in frames)
    ε = Array{ Int64 }( undef, n, 2 ); # preallocate
    ε[ :, 1 ] = collect( 1:ω:NRecFrames ); # start and
    ε[ :, 2 ] = ε[ :, 1 ] .+ ω .- 1; # end in frames of each chunk (to cut)
    if !isinteger( NRecFrames / n )
        ε[ end, 2 ] = NRecFrames; # end in frames of each chunk (to cut)
    end
    return ε, ω
end

"""
    div_ab( n::Int64, lo::Int = 1, hi::Int = n )
        Divisors of the number n between the values "lo" and "hi", if they are not defined then it takes from 1 to n
"""
function div_ab( n::Int, lo::Int = 1, hi::Int = n )
    ρ = collect( 1:floor( Int, sqrt( n ) ) ) ; # the numbers one by one, from the square root
    σ1 = findall( n.%ρ .== 0 ); # square root divisors ( remainder = 0 )
    σ2 = Int.( ( n ) ./ ( σ1 ) ); # Take out the pairs (of 100, 2-50, 10-10, etc.)
    σ = sort( unique( vcat( σ1, σ2 ) ) ); # remove duplicates, concatenate, sort
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


end

"""
    Collection of functions for analysis of BRW (HDF5) files generated with the BrainWave program from the company 3Brain.
    Laboratory 19 of the CINVESTAV in charge of Dr. Rafael Gutierrez Aguilar.
    Work developed by Isabel Romero-Maldonado (2020 - )
    isabelrm.biofisica@gmail.com
    https://github.com/LBitn
"""

module JuliaTools

using JLD
using Distributions
using BinningAnalysis
using StatsBase
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

end

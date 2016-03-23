module Resonators

using PainterQB
using PainterQB.VNA

export search_sidecoupled

function search_sidecoupled(ins::InstrumentVNA, startf::Real, stopf::Real, cutoff=0.3, span=5e5, markers=1:9)
    startf > stopf && error("Start frequency > stop frequency.")

    # Calculate frequency steps
    stepf = stopf/span
    numpts = length(startf:stepf:stopf)

    configure(ins, FrequencyStart, startf)
    configure(ins, FrequencyStop, stopf)
    configure(ins, NumPoints, 4001)
    configure(ins, IFBandwidth, 50)

    configure(ins, Graphs, [1])
    configure(ins, NumTraces, 1)
    configure(ins, VNA.LogMagnitude, 1, 1)

    shotgun(ins, markers)
    mx = search(ins,
        [MarkerSearch(:RightPeak, 1, 1, m, cutoff, VNA.Negative()) for m in markers]...)

    isnan(mx[1]) && configure(ins, Marker, 1, false)
    for m in 2:length(markers)
        (mx[m] == mx[m-1] || isnan(mx[m])) && configure(ins, Marker, m, false)
    end

    count = 0
    fs = Float64[]
    for m in markers
        if inspect(ins, Marker, m) == true
            count += 1
            push!(fs, mx[m])
        end
    end

    if count == 0
        warning("No peaks were found.")
        return
    end

    configure(ins, NumTraces, 3)
    configure(ins, Graphs, [1 3; 2 3])
    configure(ins, VNA.LogMagnitude, 1, 1)
    configure(ins, VNA.Phase, 1, 2)
    configure(ins, VNA.PolarComplex, 1, 3)
    configure(ins, FrequencyCenter, fs[1])
    configure(ins, FrequencySpan, span)

    # count == 2 && configure(ins, Graphs, [1 2])
    # 3 <= count <= 4 && configure(ins, Graphs, [1 2; 3 4])
    # 5 <= count <= 6 && configure(ins, Graphs, [1 2 3; 4 5 6])
    # 7 <= count <= 8 && configure(ins, Graphs, [1 2 3 4; 5 6 7 8])
    # count == 9 && configure(ins, Graphs, [1 2 3; 4 5 6; 7 8 9])
    # 10 <= count <= 12 && configure(ins, Graphs, [1 2 3 4; 5 6 7 8; 9 10 11 12])
    # 13 <= count && configure(ins, Graphs, [1 2 3 4; 5 6 7 8; 9 10 11 12])
    mx
end

end # module

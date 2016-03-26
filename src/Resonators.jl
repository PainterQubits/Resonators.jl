module Resonators

using PainterQB
using PainterQB.VNA
using PainterQB.E5071CModule
using DataFrames

import PainterQB.measure

export search_sidecoupled
export FrequencySweep

type FrequencySweep <: Response
    ins::InstrumentVNA
end

function measure(x::FrequencySweep)
    df = DataFrame()

    configure(x.ins, VNA.NumTraces, 3)
    configure(x.ins, VNA.Graphs, [1 2; 1 3])
    configure(x.ins, VNA.PolarComplex)
    configure(x.ins, VNA.LogMagnitude, 1, 2)
    configure(x.ins, VNA.Phase, 1, 3)
    configure(x.ins, VNA.S21, 1, 2)
    configure(x.ins, VNA.S21, 1, 3)
    configure(x.ins, VNA.S21)
    ti = inspect(x.ins, SweepTime, 1)

    configure(x.ins, BusTrigger)
    trig1(x.ins)
    sleep(ti)
    opc(x.ins)

    df[:f] = stimdata(x.ins)
    for (n,s) in [(:S11, VNA.S11), (:S12, VNA.S12), (:S21, VNA.S21), (:S22, VNA.S22)]
        configure(x.ins, s)
        df[n] = VNA.data(x.ins, VNA.PolarComplex)
    end

    df
end

function search_sidecoupled(ins::InstrumentVNA, startf::Real, stopf::Real, cutoff=3, span=5e5, markers=1:9)
    startf > stopf && error("Start frequency > stop frequency.")

    # Calculate frequency steps
    stepf = stopf/span
    numpts = length(startf:stepf:stopf)

    configure(ins, FrequencyStart, startf)
    configure(ins, FrequencyStop, stopf)
    configure(ins, NumPoints, numpts)
    configure(ins, VNA.IFBandwidth, 500)

    configure(ins, VNA.Graphs, [1])
    configure(ins, VNA.NumTraces, 1)
    configure(ins, VNA.LogMagnitude, 1, 1)
    configure(ins, VNA.S21)
    ti = inspect(ins, SweepTime)

    trig1(ins)
    sleep(ti)
    opc(ins)

    shotgun(ins, markers)
    mx = search(ins,
        [MarkerSearch(:RightPeak, 1, 1, m, cutoff, VNA.Negative()) for m in markers]...)

    isnan(mx[1]) && configure(ins, VNA.Marker, 1, false)
    for m in 2:length(markers)
        (mx[m] == mx[m-1] || isnan(mx[m])) && configure(ins, VNA.Marker, m, false)
    end

    count = 0
    fs = Float64[]
    for m in markers
        if inspect(ins, VNA.Marker, m) == true
            count += 1
            push!(fs, mx[m])
        end
    end

    if count == 0
        warning("No peaks were found.")
        return
    end

    out = DataFrame[]
    for f in fs
        configure(ins, VNA.FrequencyCenter, f)
        configure(ins, VNA.FrequencySpan, span)
        push!(out, measure(FrequencySweep(ins)))
    end

    out
    # count == 2 && configure(ins, Graphs, [1 2])
    # 3 <= count <= 4 && configure(ins, Graphs, [1 2; 3 4])
    # 5 <= count <= 6 && configure(ins, Graphs, [1 2 3; 4 5 6])
    # 7 <= count <= 8 && configure(ins, Graphs, [1 2 3 4; 5 6 7 8])
    # count == 9 && configure(ins, Graphs, [1 2 3; 4 5 6; 7 8 9])
    # 10 <= count <= 12 && configure(ins, Graphs, [1 2 3 4; 5 6 7 8; 9 10 11 12])
    # 13 <= count && configure(ins, Graphs, [1 2 3 4; 5 6 7 8; 9 10 11 12])

end

end # module

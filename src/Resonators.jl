module Resonators

using PainterQB
using PainterQB.VNA

export search_sidecoupled

function search_sidecoupled(ins::InstrumentVNA, startf::Real, stopf::Real, cutoff=0.3, resolve_q=5e5, markers=1:9)
    startf > stopf && error("Start frequency > stop frequency.")

    # Calculate frequency steps
    stepf = stopf/resolve_q
    numpts = length(startf:stepf:stopf)

    configure(ins, FrequencyStart, startf)
    configure(ins, FrequencyStop, stopf)
    configure(ins, NumPoints, 4001)
    configure(ins, IFBandwidth, 50)

    shotgun(ins, markers)
    mx = search(ins,
        [MarkerSearch(:RightPeak, 1, 1, m, cutoff, VNA.Negative()) for m in markers]...)

    for m in 2:length(mx)
        (mx[m] == mx[m-1] || isnan(mx[m])) && configure(ins, Marker, m, false)
    end




end

end # module


type Dyadic
    num::Int64
    log2den::Int64
    function Dyadic(num::Integer, den::Integer)
        tmp = log2(abs(den))
        if isinteger(tmp)
            return new( sign(den)*num, Int(tmp) )
        else
            error( DomainError() )
        end
    end
end
convert(Float64, x::Dyadic) = x.num / 2^(x.log2den)


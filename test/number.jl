#
# count the number of surreals created by generation n 
#    https://oeis.org/A174972
#


function a(n::Integer)
    a = zero(n)
    for i=1:2^n
        a += i * binomial(2^n - 1, i-1)
    end
    return a
end

# 8 is as high as OEIS goes
# approximately doubles the number of digits each time
for n=zero(BigInt) : 9
    println("  n = $n, a(n) = ", a(n), ", approx(a(n)) = ", float(a(n)))
end

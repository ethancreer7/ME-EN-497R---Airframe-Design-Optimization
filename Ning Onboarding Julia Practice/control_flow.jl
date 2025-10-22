function camber(x,p,c)
    if x <= p
        return (c*(2*p*x - x^2))/(p^2)
    else
        return (c)*(1 - 2*p + 2*p*x - x^2)/((1-p)^2)
    end
end

println("\e[4m  x | z  \e[0m")
for i in range(0.0, 1.0, step=0.1)
    z = camber(i, 0.4, 0.02)
    println("$i | $z")
end

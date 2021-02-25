module VariableTransforms

    export from_11_to_R
    export from_R_to_11
    export from_01_to_R
    export from_R_to_01
    export from_0L_to_R
    export from_R_to_0L
    export from_pos_to_R
    export from_R_to_pos


    from_11_to_R(x)=log1p(x)-log1p(-x)
    from_R_to_11(x)=2*exp(x)/(1+exp(x))-1.

    from_01_to_R(x)=log(x)-log1p(-x)
    from_R_to_01(x)=exp(x)/(1+exp(x))

    from_0L_to_R(x; fmax=1, steepness=0.1, midpoint=0)=(steepness*midpoint+log(x)-(fmax==1 ? log1p(-x) : log(fmax-x)))/steepness
    from_R_to_0L(x; fmax=1, steepness=0.1, midpoint=0)=fmax/(1+exp(-steepness*(x-midpoint)))

    from_pos_to_R(x)=log(x)
    from_R_to_pos(x)=exp(x)

    from_11_to_R_prime(x)=2/(1-x^2)
    from_01_to_R_prime(x)=1/(x-x^2)
    from_0L_to_R_prime(x; fmax=1, steepness=0.1, midpoint=0)=(1/x + 1/(fmax-x))/steepness
    from_pos_to_R_prime(x)=1/x

    from_R_to_11_prime(x)=2*from_R_to_01(x)*(1-from_R_to_01(x))
    from_R_to_01_prime(x)=from_R_to_01(x)*(1-from_R_to_01(x))
    from_R_to_pos_prime(x)=from_R_to_pos(x)
    from_R_to_0L_prime(x; fmax=1, steepness=0.1, midpoint=0)=fmax*steepness*from_R_to_0L(x; fmax=1, steepness=steepness, midpoint=midpoint)*(1-from_R_to_0L(x; fmax=1, steepness=steepness, midpoint=midpoint))
end # module


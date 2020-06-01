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

end # module

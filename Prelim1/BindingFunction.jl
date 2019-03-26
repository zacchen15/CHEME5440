function binding_function(saturation_constant,inducer_concentration,binding_order)
    K = saturation_constant
    I = inducer_concentration
    n = binding_order

    f = (I^n) / (K^n + I^n)
    return f
end

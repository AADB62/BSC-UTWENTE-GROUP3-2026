function g=g_auto(eda,eds)
    u=eda./eds;
    g=4.5./(8*u.^4+1)+0.5;
end
function gain=invTrue_rms(dBSPL)
peff=2e-5*10.^(dBSPL/20);
gain=peff/(2*10^.5);
end
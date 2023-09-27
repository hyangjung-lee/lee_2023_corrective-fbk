function vect_full = get_Fpostevent(tmp_sensIDX, tmp, tmp2)

vect_full = tmp(tmp_sensIDX); 
vect_full(tmp2(tmp_sensIDX) == 1) = vect_full(tmp2(tmp_sensIDX) == 1)*10;
vect_full(tmp2(tmp_sensIDX) == -1) = vect_full(tmp2(tmp_sensIDX) == -1)*13;

end

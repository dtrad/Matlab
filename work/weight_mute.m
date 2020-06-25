function weight_mute(it1,it2,iq1,iq2)
w=ones(512,100);
for it=it1:it2 
   for iq=iq1:iq2
        w(it,iq)=(exp(-(10-(abs(it2-it)+abs(iq2-iq)+abs(it1-it)+abs(iq1-iq)))))+1;
   end
end
figure,mesh(w)
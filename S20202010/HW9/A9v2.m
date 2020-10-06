%Assignment 9
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = 'What is the case #?';
N= input(prompt);

%%%%%%%%%%% case choice %%%%%%%%%%%%%

if (N == 1)
    body = 0;
    drain = 1;
    source = 0;
elseif (N == 2)
    body = 0;
    drain = 0;
    source = 1;
elseif (N == 3)
    body = 1;
    drain = 0;
    source = 0;
elseif (N == 4)
    body = 1;
    drain = 1;
    source = 1;
else
    disp(sprintf('wrong number'));
    
end
    

A = zeros(45,45);
b = zeros(45,1);
nx = 9 ; %number of coloum


        

%%%%%%%%%%%%% construct vector b %%%%%%%%%%%%%%%
   if (N==1|N==2|N==3|N==4)
       for iy = 1:5
            for ix = 1:9
                if (iy == 1) %substrate body 
                    b(((iy-1)*nx+ix), 1) = body;
                elseif(iy==5)
                        if(ix==1|ix==2) % source 
                            b(((iy-1)*nx+ix), 1) = source;
                        elseif(ix==8|ix==9) % drain
                            b(((iy-1)*nx+ix), 1) = drain;
                        end
                else
                end
                
                    
            end
       end
   end
   
  %%%%%%%%%%%%%%%%%%% construct A with empty circles %%%%%%%%%%%%%%%%%%%%%%
     if (N==1|N==2|N==3|N==4)
        for iy = 2:4
            for ix = 2:8
                     A(((iy-1)*nx+ix), ((iy-1)*nx+ix-1)) = 1;
                     A(((iy-1)*nx+ix), ((iy-1)*nx+ix+1)) = 1;
                     A(((iy-1)*nx+ix), ((iy-1)*nx+ix)) = -4;
                     A(((iy-1)*nx+ix), ((iy-2)*nx+ix)) = 1; 
                     A(((iy-1)*nx+ix), ((iy)*nx+ix)) = 1;   
            end
        end
     end
  
    %%%%%%%%%%%%%%%%%%% construct A for the top and bottom %%%%%%%%%%%%%%%%%%%%%%
     if (N==1|N==2|N==3|N==4)
        for iy = 1:5
            if iy == 1
                for ix = 1:9
                    A(((iy-1)*nx+ix), ((iy-1)*nx+ix)) = 1;
                end
            elseif iy == 5
                for ix = 1:9
                    A(((iy-1)*nx+ix), ((iy-1)*nx+ix)) = 1;
                end
            end
        end
     end
     
    
     
     %%%%%%%%% construct A for the Neumann boundary condition%%%%%%%%%%%%
     if (N==1|N==2|N==3|N==4)
      for iy = 2:4
             for ix = 1:9
                 if ix==1 %lefthand
                    A(((iy-1)*nx+ix), ((iy-1)*nx+ix+1)) = 1; 
                    A(((iy-1)*nx+ix), ((iy-1)*nx+ix)) = -2;
                    A(((iy-1)*nx+ix), ((iy)*nx+ix)) = 0.5;
                    A(((iy-1)*nx+ix), ((iy-2)*nx+ix)) = 0.5;
                 elseif ix==9 %right hand
                    A(((iy-1)*nx+ix), ((iy-1)*nx+ix-1)) = 1; 
                    A(((iy-1)*nx+ix), ((iy-1)*nx+ix)) = -2;
                    A(((iy-1)*nx+ix), ((iy)*nx+ix)) = 0.5;
                    A(((iy-1)*nx+ix), ((iy-2)*nx+ix)) = 0.5;
                end
        
            end
        end
     end
     
     if (N==1|N==2|N==3|N==4)
     x=A\b;
     end
     
     psi_2d=reshape(x,9,5);

     surf(psi_2d);
     xlabel('y');
     ylabel('x');
     zlabel('phi');
     
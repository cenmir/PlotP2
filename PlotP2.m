classdef PlotP2 < matlab.mixin.Copyable
    %PlotP2 Plots a P2 Triangulated surface
    %   Detailed explanation goes here
    
    properties
        light
        patch
        edge
    end
    properties (Hidden)
        X
        tri
        CData
        xnod
        ynod
        znod
    end
    
    properties (Dependent)
        EdgeColor
        FaceColor
        
    end
    
    methods
        function H = PlotP2(tri, X, np, C,varargin)
            %%
            nP2ele=size(tri,1);
            if isequal(C,0) || strcmpi(C,'none')
                CData = X(:,3);
            elseif strcmpi(C,'k')
                CData = X(:,3);
            elseif strcmpi(C,'y')
                CData = X(:,3);
            elseif strcmpi(C,'r')
                CData = X(:,3);
            elseif strcmpi(C,'g')
                CData = X(:,3);
            elseif strcmpi(C,'b')
                CData = X(:,3);
            elseif strcmpi(C,'c')
                CData = X(:,3);
            elseif strcmpi(C,'m')
                CData = X(:,3);
            elseif strcmpi(C,'w')
                CData = X(:,3);
            else
                
            end
            nTriEle = nP2ele*np^2;
            if nTriEle>30000
                warning('PLOTP2:manyelements',['Number of generates triangular elements is high (',num2str(nTriEle),')!'])
            end
            H.xnod = X(:,1); H.ynod = X(:,2); H.znod = X(:,3);
            dx=1/np;
            [X,Y] = meshgrid(0:dx:1,0:dx:1);
            x=X(:);y=Y(:);
            ind=find(y>eps+1-x);x(ind)=[];y(ind)=[];
            
            h = varargin{1};
            if nargin > 4
                xfigure(h);
            else
                xfigure;
            end
            
            axis equal;

            hold on;
            
            % Loop over all elements and map all P1 coordinates to P2 using the second
            % order or P2 triangular element base function
            [T,X,U,K,k] = H.P2elem(tri, nP2ele, np, x, y, CData);
            
            % Draw the whole mesh! With or without options
            xn=X(:,1);yn=X(:,2);zn=X(:,3);
            H.patch =  patch('faces',T,'vertices',[xn(:) yn(:) zn(:)],'CData',U(:));
            if ~isequal(C,H.CData)
                set(H.patch,'FaceColor',C)
            end
            set(H.patch,'EdgeColor','none')
            set(H.patch,'FaceLighting','gouraud')
            H.light = light;
            
            kk1=reshape(K,length(k),length(K)/length(k))';
            H.edge=trimesh(kk1,xn,yn,zn);
            set(H.edge,'FaceColor','none')
            set(H.edge,'EdgeColor','k')
            set(H.edge,'EdgeLighting','none')
            
            
        end
        
        function EdgeColor = get.EdgeColor(H)
           EdgeColor = H.edge.EdgeColor;
        end  
        function set.EdgeColor(H,C)
           H.edge.EdgeColor = C;
        end  
        
        function FaceColor = get.FaceColor(H)
           FaceColor = H.patch.FaceColor;
        end 
        function set.FaceColor(H,C)
           H.patch.FaceColor = C;
        end  
        
    end
    
    methods (Access = private)
        
        function [T,X,U,K,k] = P2elem(H,tri, nele, np, x, y, C)
            % The following is done to go from plotting one P2 element at a time
            % to all of them at once. H is done by putting all local tt into T, all
            % local coordinates [xn,yn,zn] into X and all uc3 into U
            T = zeros(nele*np^2,3);
            X = zeros(nele*length(x)^2,3);
            U = zeros(nele*length(x)^2,1);
            for iel=1:nele
                iv6=tri(iel,:);
                iv3=tri(iel,1:3);
                xc=H.xnod(iv6);yc=H.ynod(iv6);zc=H.znod(iv6);
                
                
                uc3=C(iv3);
                if size(C,1)==3 %usol can be sent in as [data,data,data]'
                    uc3=C(:,iel);
                end
                
                
                fiP2=[1-3*x+2*x.^2-3*y+4*x.*y+2*y.^2,x.*(-1+2*x),y.*(-1+2*y),-4*x.*(-1+x+y),4*x.*y,-4*y.*(-1+x+y)];
                fiP1=[1-x-y,x,y];
                
                xn=fiP2*xc;
                yn=fiP2*yc;
                zn=fiP2*zc;
                
                % Insert coordinates from iel into X, that will contain all coordinates
                lx=length(xn);
                if iel==1
                    X(iel:lx,:)=[xn,yn,zn];
                else
                    X((iel-1)*lx+1:iel*lx,:)=[xn,yn,zn];
                end
                
                %     unP2=fiP2*uc6;% For later, when solutions needs to be P2
                
                unP1=fiP1*uc3;
                % Global solution vector is assembled
                lu=length(unP1);
                if iel==1
                    U(iel:lu,1)=unP1;
                else
                    U((iel-1)*lu+1:iel*lu,:)=unP1;
                end
                
                
                % replaced by subTriang, which is faster but specialised for H purpose
                
                tt = H.subTriang(np,x);
                
                
                % Global triangle matrix is assembled
                lt=size(tt,1);
                if iel==1
                    T(iel:lt,:)=tt;
                    tt0=tt;
                    tt1=tt0;
                else
                    maxtt=max(tt(:));
                    tt1=tt1+maxtt;
                    lt=size(tt1,1);
                    T((iel-1)*lt+1:iel*lt,:)=tt1;
                end
                
                % Global element edge vector is assembled
                % H is needed to draw the P2 element edges
                indx=find(x == min(x));
                indy=find(y == min(y));
                ns=np+1;
                ns3=zeros(1,ns);
                ns3(1,1)=ns;
                ii=2;
                for i1=ns-1:-1:1
                    ns3(1,ii)=ns3(1,ii-1)+i1;
                    ii=ii+1;
                end
                k = [indx;ns3(2:end)';indy(end-1:-1:1)];
                
                lk = length(k);
                if iel==1
                    kk(iel:lk,1)=k;
                    k0=k;
                    k1=k0;
                else
                    k2=k1+max(k);
                    k0=k1;
                    k1=k2;
                    kk((iel-1)*lk+1:iel*lk,1)=k1;
                end
                K = kk;
                
            end
            
        end
        
        function tt = subTriang(H,n,x)
            nT=length(x);
            ns=n+1;
            
            
            ns3=zeros(1,ns);
            ns3(1,1)=ns;
            ii=2;
            for i1=ns-1:-1:1
                ns3(1,ii)=ns3(1,ii-1)+i1;
                ii=ii+1;
            end
            % ns3; %Third side T
            
            emax=-1;
            for i1=1:n
                emax=emax+2;
            end
            % emax; %number T in last sector
            elements=1:2:emax;
            nelem = sum(elements); %number elements
            T1 = 1:nT; %number T
            
            si=1;
            ei=ns;
            T2=zeros(ns);
            for i1=1:ns
                T2(i1,i1:ns) = T1(si:ei);
                si=ei+1;
                ei=ns3(i1)+ns-i1;
            end
            
            tri2=zeros(nelem,3);
            i1=1;
            giel=1;
            for iseq=elements
                if iseq==1
                    seq1 = T2(:,[1,2]);
                    ind1= seq1==0;
                    seq1(ind1)=[];
                    seq1 = [seq1(3),seq1(2),seq1(1)];
                    tri2(1,:)=seq1;
                    giel=2;
                end
                if iseq~=1
                    seqi = T2(:,[i1,i1+1]);
                    seqi1=seqi(:,1);
                    seqi2=seqi(:,2);
                    ind1= seqi1==0;
                    seqi1(ind1)=[];
                    c1=length(seqi1);
                    
                    ind2= seqi2==0;
                    seqi2(ind2)=[];
                    c2=length(seqi2);
                    
                    seqi=seqi(:);
                    ind1= seqi==0;
                    seqi(ind1)=[];
                    seqi = seqi';
                    a = c2-1;
                    b = iseq - a;
                    for ai=1:a
                        tri2(giel,:) = seqi([ai,ai+c1+1,ai+c1]);
                        giel=giel+1;
                    end
                    for bi=1:b
                        tri2(giel,:) = seqi([bi,bi+1,bi+c1+1]);
                        giel=giel+1;
                    end
                end
                i1=i1+1;
            end
            tt=tri2;
            
            
            
        end
        
        
        
    end
    
end


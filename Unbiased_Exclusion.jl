using Plots, SpecialFunctions, Random
gr()

function Stochastic(LX,LY,Tplot,t1,A0,MC,PM)
    Q=Int(sum(A0))
    T=Int(t1)
    AA=zeros(LX,LY)
    density=zeros(LX)
    pos0=zeros(Q,2)
    post=zeros(Q,2)
    
    
    agent=0
    for i in 1:LX
        for j in 1:LY
            if A0[i,j] == 1.0
            agent = agent+1
            pos0[agent,1]=i-LX/2
            pos0[agent,2]=j
            end 
        end
    end
    
    for MM in 1:MC
        println(MM)
        AAtemp=copy(A0)
        for kk in 1:T
        
    
            count = 0
            while count  < Q
                II =rand(1:LX)
                JJ=rand(1:LY)
                if AAtemp[II,JJ] > 0.0 && II>1 && II<LX && JJ > 1 && JJ < LY 
                count=count+1
                R =rand(1)
                S=rand(1)
                    if AAtemp[II-1,JJ] == 0.0 && R[1] > 0 && R[1]<=1/4 && S[1] <=PM
                    AAtemp[II,JJ]=0.0
                    AAtemp[II-1,JJ]=1.0
                    elseif AAtemp[II+1,JJ] == 0.0 && R[1] > 1/4 && R[1]<=2/4  && S[1] <=PM
                    AAtemp[II,JJ]=0.0
                    AAtemp[II+1,JJ]=1.0
                    elseif AAtemp[II,JJ-1] == 0.0 && R[1] > 2/4 && R[1]<=3/4  && S[1] <=PM
                    AAtemp[II,JJ]=0.0
                    AAtemp[II,JJ-1]=1.0
                    elseif AAtemp[II,JJ+1] == 0.0 && R[1] > 3/4 && R[1]<=4/4  && S[1] <=PM
                    AAtemp[II,JJ]=0.0
                    AAtemp[II,JJ+1]=1.0
                    end
    
                    elseif AAtemp[II,JJ] > 0.0 && II>1 && II<LX && JJ == 1 
                        count=count+1
                        R =rand(1)
                        S=rand(1)
                            if AAtemp[II-1,JJ] == 0.0 && R[1] > 0 && R[1]<=1/4 && S[1] <=PM
                            AAtemp[II,JJ]=0.0
                            AAtemp[II-1,JJ]=1.0
                            elseif AAtemp[II+1,JJ] == 0.0 && R[1] > 1/4 && R[1]<=2/4 && S[1] <=PM
                            AAtemp[II,JJ]=0.0
                            AAtemp[II+1,JJ]=1.0
                            elseif AAtemp[II,JJ+1] == 0.0 && R[1] > 3/4 && R[1]<=4/4 && S[1] <=PM
                            AAtemp[II,JJ]=0.0
                            AAtemp[II,JJ+1]=1.0
                            end
            
    
                        elseif AAtemp[II,JJ] > 0.0 && II>1 && II<LX &&  JJ == LY 
                                count=count+1
                                R =rand(1)
                                S=rand(1)
                                    if AAtemp[II-1,JJ] == 0.0 && R[1] > 0 && R[1]<=1/4  && S[1] <=PM
                                    AAtemp[II,JJ]=0.0
                                    AAtemp[II-1,JJ]=1.0
                                    elseif AAtemp[II+1,JJ] == 0.0 && R[1] > 1/4 && R[1]<=2/4 && S[1] <=PM
                                    AAtemp[II,JJ]=0.0
                                    AAtemp[II+1,JJ]=1.0
                                    elseif AAtemp[II,JJ-1] == 0.0 && R[1] > 2/4 && R[1]<=3/4 && S[1] <=PM
                                    AAtemp[II,JJ]=0.0
                                    AAtemp[II,JJ-1]=1.0
                                     end
                           
                    
    
    
    
                end
            end
    
    
        
        if kk == Int(Tplot[2]) && MM == 1
        agent=0
        for i in 1:LX
            for j in 1:LY
                if AAtemp[i,j] > 0.0 
                agent = agent+1
                post[agent,1]=i-LX/2
                post[agent,2]=j
                end 
            end
        end    
    end
    
        end
    AA=AA+AAtemp
    
    
    end
    
    AA=AA/MC
    
    
    
    for i in 1:LX
     density[i]=0
        for j in 1:LY
    density[i]=density[i]+AA[i,j]
     end
    end
    
    density=density/LY
    
    
    return pos0,post,density
    
    
    
    end



 

LX=400
LY=50
PM=1.0
D=PM/4
T=1000
Tplot=[0,T]
MC=10
h=50

A0=zeros(LX,LY)
xxloc=zeros(LX)
yyloc=zeros(LY)
for i in 1:LX
    xxloc[i]=-LX/2+(i-1)
    for j in 1:LY
    yyloc[j]=0+(j-1)
        if abs(xxloc[i]) <= h
        A0[i,j]=1.0
        end
    end
end


C(x)=0.5*(erf((h-x)/sqrt(4*D*T))+erf((h+x)/sqrt(4*D*T)));
(pos0,pos1,CS)=Stochastic(LX,LY,Tplot,T,A0,MC,PM);



p1=scatter(xxloc,CS,mc=:blue,msc=:match,label=false)
p1=plot!(C,-LX/2,LX/2,lw=3,lc=:red,label=false)
q1=plot!(xlims=(-LX/2,LX/2),xticks=([-200,-100,0,100,200],[L"-200",L"-100",L"0",L"100",L"200"]))
q1=plot!(ylims=(0,1.00),yticks=([0,0.25,0.5,0.75,1.0],[L"0.00", L"0.25", L"0.50", L"0.75", L"1.00"]))
q1=plot!(xguidefontsize=12, yguidefontsize=12,xtickfontsize=12, ytickfontsize=12)
q1=plot!(xlabel=L"x",ylabel=L"c(x,t)")
display(p1)


q0=scatter(pos0[:,1],pos0[:,2],markersize=2,markershape=:circle, markercolor=:blue,msw=0,legend=false,aspect_ratio=:equal)
q0=plot!(xlims=(-LX/2,LX/2),xticks=([-200,-100,0,100,200],[L"-200",L"-100",L"0",L"100",L"200"]))
q0=plot!(ylims=(0,h),yticks=([0,h],[L"0",L"50"]))
q1=scatter(pos1[:,1],pos1[:,2],markersize=2,markershape=:circle, markercolor=:blue,msw=0,legend=false,aspect_ratio=:equal)
q1=plot!(xlims=(-LX/2,LX/2),xticks=([-200,-100,0,100,200],[L"-200",L"-100",L"0",L"100",L"200"]))
q1=plot!(ylims=(0,h),yticks=([0,h],[L"0",L"50"]))
q2=plot(q0,q1,layout=(2,1))
display(q2)   


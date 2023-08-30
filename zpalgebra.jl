#-------------------------------------------------
#-------------------------------------------------

function rank_p!(A,p)
    # find rank(A) in Z_p
    # changes A 
    m,n=size(A)
    lastpivotcol=0

    for i in 1:m
        for j in (lastpivotcol+1):n
            if A[i,j]!=0
                for jj in (j+1):n
                    if A[i,jj]!=0
                        A[i:end,jj]=@views mod.( A[i:end,jj].+  A[i:end,j].*(multiplicativeInverse(A[i,j],p)*(p-A[i,jj])),p)
                    end
                end
                swapcol!(A,j,lastpivotcol+1,i)
                lastpivotcol+=1
                break
            end
        end
    end

    return lastpivotcol
end
#-------------------------------------------------
#-------------------------------------------------
function swapcol!(A,i,j,toprow=1)
    #swap columns i and j
    #if toprow is given it only swaps the A[toprow:end,i] with A[toprow:end,j]
    # doesn't depend on Zp algebra

    if i!=j
        temp=A[toprow:end,i]
        A[toprow:end,i]=@view A[toprow:end,j]
        A[toprow:end,j]=temp
    end
    return
end
#-------------------------------------------------
#-------------------------------------------------
function swaprow!(A,i,j,leftcol=1)
    #swap rows i and j
    #if toprow is given it only swaps the A[i,leftcol:end] with A[j,leftcol:end]
    # doesn't depend on Zp algebra

    if i!=j
        temp=A[i,leftcol:end]
        A[i,leftcol:end]=@view A[j,leftcol:end]
        A[j,leftcol:end]=temp
    end
    return
end

#-------------------------------------------------
#-------------------------------------------------


function extendedEuclideanAlg(a,b,p) 
    # Extended Euclidean algorithm to solve for r=gcd(a,b),s,t 
    # in a*s+b*t=r for given a and b, assuming the filed to be Z_p 
    # for some prime(?) p.

        T=typeof(a)
    
    if Int(p)^2>typemax(T)
        println("OVERFLOW WARNING: p^2 > typemax")
    end


    r0=a
    r1=b

    s0=T(1)
    s1=T(0)

    t0=T(0)
    t1=T(1)


    r=r0%r1



    q=(r0-r)÷r1
    
    s=mod(s0-q*s1,p)
    t=mod(t0-q*t1,p)



    while r != 0
        temp=r1
        r1=r
        r0=temp

        r=r0%r1
        q=(r0-r)÷r1

        s0=s1
        s1=s
        
        t0=t1
        t1=t
        
        s=mod(s0-q*s1,p)
        t=mod(t0-q*t1,p)
        
    end

    return r1,s1,t1
end

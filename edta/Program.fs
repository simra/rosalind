﻿// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.
open System
open System.IO

let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()
let splitNewline (x:string) = x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)

type FASTA =
    {Label:string;String:string}
    
let parseFasta (text:string) =
    text.Split([|'>'|],StringSplitOptions.RemoveEmptyEntries)
    |> Seq.map 
        (fun lines-> 
            lines.Split([|'\n';'\r'|],StringSplitOptions.RemoveEmptyEntries)
            |> fun toks -> 
                {Label=toks.[0]; String=toks.[1..]|>String.concat ""})

let guard m c =
    if Map.containsKey c m then
        m.[c]
    else
        0

let rec partFact (n:int64) (stop:int64) =
    if (n<=stop) then 1L
    else n*(partFact (n-1L) stop)

let rec fact (n:int64) = partFact n 1L 

let rec partFactB (n:bigint) (stop:bigint) =
    if (n<=stop) then bigint 1
    else n*(partFactB (n-bigint 1) stop)

let rec factB (n:bigint) = partFactB n (bigint 1) 

// https://en.wikipedia.org/wiki/Binomial_coefficient#Multiplicative_formula
// This was actually incorrect on the wikipedia page and I had to make an edit!
let rec C n k =
    if k=0L || k=n then 1L
    else if k>=(n/2L+1L) then C n (n-k)
    else 
        [1L..k]
        |> Seq.map (fun i -> (float (n-k+i))/(float i))
        |> Seq.reduce (*)
        |> int64

open System.Collections.Generic
 // some handy dynamic programming snippets from http://www.fssnip.net/8P
 /// The function creates a function that calls the argument 'f'
 /// only once and stores the result in a mutable dictionary (cache)
 /// Repeated calls to the resulting function return cached values.
let memoize f =    
   // Create (mutable) cache that is used for storing results of 
   // for function arguments that were already calculated.
    let cache = new Dictionary<_, _>()
    (fun x ->
       // The returned function first performs a cache lookup
       let succ, v = cache.TryGetValue(x)
       if succ then v else 
         // If value was not found, calculate & cache it
         let v = f(x) 
         cache.Add(x, v)
         v)



// EDTA
let (^^) a b = (a && not b) || (not a && b)
let aligned (s:string) (t:string) =
    s.Length=t.Length &&
        Seq.zip s t 
        |> Seq.map (fun (a,b) -> a=b || ((a='-') ^^ (b='-') ))
        |> Seq.reduce (&&)

let edta (s:string) (t:string) = 
    let rpt x = [for i in 1..x do yield "-"]|>String.concat ""   
    let rec helper = memoize (fun (i:int,j:int) ->
        //eprintfn "%d %d" i j
        if i = 0 then
            if j = 0 then (0,"","")
            else (j,rpt j,t.Substring(0,j))
        else if j = 0 then
            (i,s.Substring(0,i),rpt i)
        else    
            seq {
                let (a1,b1,c1)= helper(i-1,j)
                yield (1+a1,b1+string s.[i-1],c1+"-")
                let (a2,b2,c2)= helper(i,j-1)
                yield (1+a2,b2+"-",c2+string t.[j-1])
                let (a3,b3,c3)= helper(i-1,j-1)
                let di=if s.[i-1]=t.[j-1] then 0 else 1
                yield (di+a3,b3+string s.[i-1],c3+string t.[j-1])
            } 
            |> Seq.minBy (fun (a,b,c)->a))
    helper(s.Length,t.Length)

    (*
    if aligned s t then (0,s,t)
    else if s="" then 
        let (a,b,c)= edta("-",t)
        (1+a,b,c)
    else if t="" then 
        let (a,b,c) =edta(s,"-")
        (1+a,b,c)
    else 
        seq {
            if not (s.[0]='-' && t.[0]='-') then
                let di=if s.[0]=t.[0] then 0 else 1
                let (a,b,c)=edta(s.Substring(1),t.Substring(1))
                yield (di+a,string s.[0]+b,string t.[0]+c)
               (* if (s.[0]<>'-' && s.Length<t.Length) then
                    let (a,b,c) = edta("-"+s,t)
                    yield (1+a,b,c)
                if (t.[0]<>'-' && t.Length<s.Length) then 
                    let (a,b,c) = edta(s,"-"+t)
                    yield (1+a,b,c) *)            
            else 
                if t.[0]<>'-' then
                    let (a,b,c)=edta("-"+s,t)
                    yield (1+a,b,c)
                if s.[0]<>'-' then
                    let (a,b,c)=edta(s,"-"+t)
                    yield (1+a,b,c)
          }
        |> Seq.minBy (fun (a,b,c)->a)
        )
        *)



let ctea (s:string) (t:string) = 
    let rpt x = [for i in 1..x do yield "-"]|>String.concat ""   
    let rec helper = memoize (fun (i:int,j:int) ->
        //eprintfn "%d %d" i j
        if i = 0 then
            if j = 0 then (0,1L)
            else (j,1L)
        else if j = 0 then
            (i,1L)
        else    
            seq {
                let (l1,c1)= helper(i-1,j)
                yield (1+l1,c1)
                let (l2,c2)= helper(i,j-1)
                yield (1+l2,c2)
                let (l3,c3)= helper(i-1,j-1)
                let di=if s.[i-1]=t.[j-1] then 0 else 1
                yield (di+l3,c3)
            } 
            |> Seq.groupBy (fun (l,c)->l)
            |> Seq.minBy (fun (l,s)->l)
            |> fun (l,s) -> (l,s|>Seq.fold (fun ci (l,c) -> (ci+c)%134217727L) 0L)
            )
    helper(s.Length,t.Length)



let BLOSUM62 =
    @"A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7"
    |> splitNewline
    |> fun toks -> 
        let alph=
            toks.[0].Split([|' '|],StringSplitOptions.RemoveEmptyEntries)
            |> Seq.mapi (fun i c -> (i+1,c.[0]))
            |> Map.ofSeq
        toks.[1..toks.Length-1]
        |> Seq.fold (fun (m:Map<char*char,int>) (l:string) ->
            let arr=l.Split([|' '|],StringSplitOptions.RemoveEmptyEntries)
            let c1=arr.[0].[0]
            arr.[1..arr.Length-1]
            |> Seq.mapi (fun i v -> (i+1,int v))
            |> Seq.fold (fun m' (j,v) -> Map.add (c1,alph.[j]) v m') m) Map.empty
        

let glob (S:Map<char*char,int>) (s:string) (t:string) = 
    let g = -5
    let rpt x = [for i in 1..x do yield "-"]|>String.concat ""   
    let rec helper = memoize (fun (i:int,j:int) ->
        //eprintfn "%d %d" i j
        if i = 0 then
            if j = 0 then (0,"","")
            else (j*g,rpt j,t.Substring(0,j))
        else if j = 0 then
            (i*g,s.Substring(0,i),rpt i)
        else    
            seq {
                let (a1,b1,c1)= helper(i-1,j)
                yield (g+a1,b1+string s.[i-1],c1+"-")
                let (a2,b2,c2)= helper(i,j-1)
                yield (g+a2,b2+"-",c2+string t.[j-1])
                let (a3,b3,c3)= helper(i-1,j-1)
                let di=S.[(s.[i-1],t.[j-1])]
                yield (di+a3,b3+string s.[i-1],c3+string t.[j-1])
            } 
            |> Seq.maxBy (fun (a,b,c)->a))
    helper(s.Length,t.Length)

let doEdta()=
    getData "edta"
    |> parseFasta
    |> Array.ofSeq
    |> fun toks -> 
        let s=toks.[0].String
        let t=toks.[1].String
        edta s t 
    |> fun (a,b,c) -> 
        printfn "%d" a
        printfn "%s" b
        printfn "%s" c                 
        
let doCtea() =      
    getData "ctea"
    |> parseFasta
    |> Array.ofSeq
    |> fun toks -> 
        let s=toks.[0].String
        let t=toks.[1].String
        ctea s t 
    |> fun (l,c)->printfn "%d" c 

let doGlob() =
    getData "glob"
    |> parseFasta
    |> Array.ofSeq
    |> fun a -> glob BLOSUM62 a.[0].String a.[1].String
    |> fun (cost,_,_) -> printfn "%d" cost


[<EntryPoint>]
let main argv = 
    doGlob()

    0 // return an integer exit code

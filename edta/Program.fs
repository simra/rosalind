// Learn more about F# at http://fsharp.org
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


[<EntryPoint>]
let main argv = 
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
    0 // return an integer exit code

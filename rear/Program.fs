// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.
open System
open System.IO
open System.Collections.Generic


let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()
let splitNewline (x:string) = x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)

// let's just do bread-and-butter bfs
// tree structure to track probed permutations. A little sloppy but I think this will work.
// to try next: work with bigints

(*
type State =
    | Empty
    | Elements of Dictionary<int,State>
    with 
    static member add s (a:int[]) (i:int) =
        if i=a.Length then s
        else 
            match s with 
            | Elements(r) ->
                if not (r.ContainsKey(a.[i])) then
                    r.[a.[i]]<-Elements(new Dictionary<int,State>())
                State.add r.[a.[i]] a (i+1)
            | Empty -> raise (new Exception "bad add")       

    static member contains s (a:int[]) (i:int) =
        if i=a.Length then true
        else
            match s with
            | Empty -> false
            | Elements(r) ->
                if not (r.ContainsKey(a.[i])) then false
                else 
                    State.contains r.[a.[i]] a (i+1)
                    *)
let toBigint (a:int[]) =
    let mutable result = bigint 0
    for x in a do
        result<-result*(bigint 10)+bigint x
    result

// in the end int64 would have been fine.
type State = HashSet<bigint>


let reversal (a:int[]) (i:int) (j:int) =
    let result=Array.copy a
    for x in [i..j-1] do
        result.[x]<-a.[j-1-(x-i)]
    result

let differs (a:int[]) (b:int[]) =
    let mutable i=0
    while i<a.Length && a.[i]=b.[i] do        
        i<-i+1
    i<a.Length

let rear4 (target:int[]) (s:int[]) = 
    //let mutable state=Elements(new Dictionary<int,State>())//Set.add (s|>toString) Set.empty
    let mutable state = new State()
    let mutable iter=0
    let allperms =
        seq {
            for i in [0..Array.length s-2] do
                for j in [i+2..Array.length s] do
                    yield (i,j)
        }
        |> List.ofSeq
    let lex =
        let s=Array.create target.Length 0
        for i in [0..s.Length-1] do
            s.[target.[i]-1]<-i
        s

    let rec bfsLoop (q0:Queue<int*int[]>) =
        iter<- iter+1
        if q0.Count=0 then
            raise (new Exception "No solution")
        else 
            let (d,s') = q0.Dequeue()
            if (iter%10000)=0 then
                eprintfn "%d %d %d" iter d q0.Count
        
            if differs s' target |> not then d
            else   
                seq {
                    for (i,j) in allperms do
                        if lex.[s'.[i]-1]>lex.[s'.[j-1]-1] then
                            let r=reversal s' i j
                            let br=toBigint r
                            if (not (state.Contains(br))) then
                                state.Add(br)|>ignore
                                yield r    
                }
                 
                |> Seq.iter ( // still needs short circuit.
                    fun r -> q0.Enqueue (d+1,r) )
                
                (*
                allperms // not clear this speeds things up..
                |> Seq.map 
                    (fun (i,j) ->
                        if lex.[s'.[i]-1]>lex.[s'.[j-1]-1] then // only expand things that swap out-of-order values.
                            Some (async { 
                                let r = reversal s' i j
                                return (r,state.Contains(toBigint r))//State.contains state r 0)
                            })
                        else None)
                |> Seq.choose id
                |> Async.Parallel
                |> Async.RunSynchronously
                |> Seq.filter (fun (_,b)->not b)
                |> Seq.iter (fun (r,_) -> 
                    //State.add state r 0 |> ignore
                    state.Add(toBigint r)|>ignore
                    // how to short-circuit out?
                    q0.Enqueue (d+1,r) ) *)
                bfsLoop q0
    let q=new Queue<int*int[]>();
    q.Enqueue(0,s)
    bfsLoop q


[<EntryPoint>]
let main argv = 
    getData "rear"
    |> splitNewline
    |> Seq.map (fun s -> s.Split ' '|> Seq.map int|>Array.ofSeq)
    //|> Seq.take 2
    |> Array.ofSeq
    |> fun a ->
        seq {
            for i in [0 .. 2 .. a.Length-1] do
                yield async { return rear4 a.[i] a.[i+1] }
            }
    |> Async.Parallel
    |> Async.RunSynchronously
    |> Seq.map string
    |> String.concat " "
    |> printfn "%s"
    0 // return an integer exit code

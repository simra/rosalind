// Learn more about F# at http://fsharp.net
// See the 'F# Tutorial' project for more help.
open System
open System.IO

let guard m c =
    if Map.containsKey c m then
        m.[c]
    else
        0

// 2. Transcribing DNA into RNA
let dnaToRna (x:string) : string = x.Replace('T','U')

// 3. reverse complement
let complement c =
    match c with
    | 'A' -> 'T'
    | 'T' -> 'A'
    | 'C' -> 'G'
    | 'G' -> 'C'
    | _ -> raise (new Exception("Invalid nucleotide"))

let reverseComplement x =
    x
    |> Seq.map complement 
    |> List.ofSeq
    |> List.rev
    |> String.Concat

// 4. fibonnaci rabbits
let rec pairsAfterN n k =
    if n=1 || n=2 then 1L
    else 
        let current=pairsAfterN (n-1) k
        let offspring=k*(pairsAfterN (n-2) k)
        current+offspring

// 5. mortal fib rabbits
// todo: handle M
let rec pairsAfterNwithM n m =
    if n=1 then [|1L;0L;0L|]
    else if n=2 then [|0L;1L;0L|]
    else
        let lastMonth=
            pairsAfterNwithM (n-1) m
        [|lastMonth.[1]+lastMonth.[2];lastMonth.[0];lastMonth.[1]|]

[<EntryPoint>]
let main argv = 
    File.ReadAllText(argv.[0])
    |> fun x -> x.Trim()
    |> Seq.countBy (fun x->x)
    |> Map.ofSeq
    |> fun m -> printfn "%d %d %d %d" (guard m 'A') (guard m 'C') (guard m 'G') (guard m 'T')
    0 // return an integer exit code





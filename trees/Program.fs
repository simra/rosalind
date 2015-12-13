// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.
open System
open System.Collections.Generic
open System.IO

let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()
let splitNewline (x:string) = x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)


type NwckTree =
    | Leaf of string
    | Anon
    | Internal of (string*NwckTree list)

let ex = fun s -> new Exception(s)

let rec parseTree (str:string) =
    if str="" then raise (ex "Empty string")
    else if str.[0]<>'(' then raise (ex "Expected open brace")
    else
        let findClose () =
            let mutable o=0
            let mutable i=1
            while (i<str.Length && (str.[i]<>')' || o>0)) do
                if str.[i]='(' then o<- o+1
                else if str.[i]=')' then o<-o-1
                i<-i+1
            i
        let findLbl i = // TODO: I think we should never see close brace here.
            if i+1>=str.Length || str.[i+1]=';' || str.[i+1]=')' || str.[i+1]=',' then ("",i+1)
            else
                let mutable j=i+2
                while j<str.Length && str.[j]<>';' && str.[j]<>',' && str.[j]<>')' do j<-j+1
                let start=i+1
                let len=j-i-1
                //eprintfn "%s %d %d %d" str str.Length start len                
                (str.Substring(start,len),j+1)
        let i = findClose ()
        let (l,j) = findLbl i
        Internal(l,parseList (str.Substring(1,i-1)))// :: parseTree (str.Substring(j))
and parseList str =
    let commaSep ()=
        Seq.unfold(
            fun (o,i,last) ->
                if i>str.Length then None
                else if i=str.Length then
                    if o=0 then
                        Some(Some(str.Substring(last+1)),(0,i+1,i))
                    else
                        raise (ex "Unexpected end of string")
                else if o=0 && str.[i]=',' then
                    Some(Some(str.Substring(last+1,i-last-1)),(0,i+1,i))
                else if str.[i]='(' then 
                    Some(None,(o+1,i+1,last))
                else if str.[i]=')' then
                    Some(None,(o-1,i+1,last))
                else Some(None,(o,i+1,last))
            ) (0,0,-1)
        |> Seq.choose id
    commaSep() 
    |> Seq.map (fun s ->
        if s="" then Anon
        else if s.[0]='(' then parseTree s
        else Leaf(s)
    )
    |> List.ofSeq

let adjacencies _nwk =
    let d = new Dictionary<string,string list>()
    let rec adj nwk : unit =        
        let checkAnon s =
            if s="" then 
                let lbl=sprintf "anon_%d" (d.Count+1)
                d.[lbl]<-List.empty
                lbl
            else 
                if not (d.ContainsKey(s)) then d.[s]<-List.empty
                s

        match nwk with
        | Leaf(s) -> ()
        | Anon -> ()
        | Internal(l,n) ->
            let l' = checkAnon l                
            n
            |> Seq.iter 
                (fun nn ->
                    match nn with
                    | Leaf(s) -> 
                        let s'=checkAnon s
                        d.[s']<-l'::d.[s']
                        d.[l']<-s'::d.[l']
                    | Anon ->
                        let s'= checkAnon ""
                        d.[s']<-l'::d.[s']
                        d.[l']<-s'::d.[l']
                    | Internal(s,sub) ->
                        let s'=checkAnon s
                        d.[s']<-l'::d.[s']
                        d.[l']<-s'::d.[l']
                        adj (Internal(s',sub))
                
                )   
    adj _nwk
    d
// inefficient?  Proper handling of dead ends?
let distance (adj:Dictionary<string,string list>) a b =
    let visited=new HashSet<string>()
    let safeValue=adj.Count*10;
    let rec dist x =
        //eprintfn "Visiting %s" x        
        visited.Add(x)|> ignore        
        if x=b then 0
        else 
           // screams out for a computation expression     
           adj.[x]
           |> Seq.map                 
                (fun n ->
                    if not (visited.Contains(n)) then
              //          eprintfn "Trying %s" n
                        1 + (dist n)
                    else
                        safeValue)
           |> Seq.min 
                     
    dist a    

getData "nwck"
|> splitNewline
|> Array.ofSeq
|> fun s -> 
        seq {
            for i in [0..2..s.Length-1] do
                yield (parseTree s.[i]|>adjacencies,s.[i+1].Split(' '))
            }
//|> Seq.take 1
|> Seq.map (fun (adj,toks) -> distance adj toks.[0] toks.[1])
|> Seq.map string
|> String.concat " "
|> printfn "%s"    

[<EntryPoint>]
let main argv = 
    printfn "%A" argv
    0 // return an integer exit code

// Learn more about F# at http://fsharp.net
// See the 'F# Tutorial' project for more help.
open System
open System.IO

let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()

// 1. DNA
let guard m c =
    if Map.containsKey c m then
        m.[c]
    else
        0

getData "dna"
|> Seq.countBy (fun x->x)
|> Map.ofSeq
|> fun m -> printfn "%d %d %d %d" (guard m 'A') (guard m 'C') (guard m 'G') (guard m 'T')

// 2. RNA: Transcribing DNA into RNA
let dnaToRna (x:string) : string = x.Replace('T','U')

getData "rna"
|> dnaToRna
|> printfn "%s"


// 3. REVC: reverse complement
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

getData "revc"
|> fun x->x.Trim()
|> reverseComplement
|> printfn "%s"

// 4. FIB: fibonnaci rabbits
let rec pairsAfterN n k =
    if n=1 || n=2 then 1L
    else 
        let current=pairsAfterN (n-1) k
        let offspring=k*(pairsAfterN (n-2) k)
        current+offspring

getData "fib"
|> fun l -> l.Split(' ')
|> fun toks -> pairsAfterN (Int32.Parse(toks.[0])) (Int64.Parse(toks.[1]))
|> printfn "%d"


// 5. FIBD: mortal fib rabbits
// todo: handle M
let rec pairsAfterNwithM n m =
    if n=1 then [1L]
    else if n=2 then [0L;1L]
    else
        let lastMonth=
            pairsAfterNwithM (n-1) m            
        (lastMonth|>Seq.skip 1|>Seq.sum) :: (lastMonth|>Seq.truncate (m-1)|>List.ofSeq)

getData "fibd"
|> fun l -> l.Split(' ')
|> fun toks -> pairsAfterNwithM (Int32.Parse(toks.[0])) (Int32.Parse(toks.[1]))
|> Seq.sum
|> printfn "%d"


// 6. GC: computing GC content
type FASTA =
    {Label:string;String:string}
    
let parseFasta (text:string) =
    text.Split([|'>'|],StringSplitOptions.RemoveEmptyEntries)
    |> Seq.map 
        (fun lines-> 
            lines.Split([|'\n';'\r'|],StringSplitOptions.RemoveEmptyEntries)
            |> fun toks -> 
                {Label=toks.[0]; String=toks.[1..]|>String.concat ""})

let gcContent s =
    let m = Seq.countBy (fun x -> x) s |> Map.ofSeq
    100.*(((guard m 'C')+(guard m 'G'))|>float)/(Seq.length s|>float)


getData "gc"
|> fun x -> x.Trim()
|> parseFasta
|> Seq.maxBy (fun f -> gcContent f.String)
|> fun f -> printfn "%s\n%f" f.Label (gcContent f.String)


// 7. HAMM: counting point mutations (hamming distance)
let dh s1 s2 =
    Seq.zip s1 s2
    |> Seq.map (fun (c1,c2) -> if c1=c2 then 0 else 1)
    |> Seq.sum


getData "hamm"
|> fun x -> x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)
|> fun lines -> printfn "%d" (dh lines.[0] lines.[1])

// 8. IPRB: Mendel's first law.
let rec partFact (n:int64) (stop:int64) =
    if (n<=stop) then 1L
    else n*(partFact (n-1L) stop)

let rec fact (n:int64) = partFact n 1L 
let C n k = (partFact n (n-k))/(fact k) // todo: check n-k>k

let pDom k m n =
    let kk = float (C k 2L)
    let km = float (k*m)
    let mm= 0.75*float (C m 2L)
    let kn= float (k*n)
    let mn=0.5*float (m*n)
    (kk+km+mm+kn+mn)/float (C (k+m+n) 2L)
// k: num  homozygous dominant
// m: num heterozygous
// n: num homozygous recessive

let p64 x = Int64.Parse(x)

getData "iprb"
|> fun l -> l.Split(' ')
|> fun t -> pDom (p64 t.[0]) (p64 t.[1]) (p64 t.[2])
|> printfn "%f"

// 9. PROT Translating RNA into Protein
// be sure source is materialized or you'll be re-traversing the sequence over and over...
let chop segmentSize source =
    source    
    |> Seq.unfold(fun chopSource ->
        if Seq.isEmpty chopSource then None else
        let segment = chopSource |> Seq.truncate segmentSize
        let rest = chopSource |> Seq.skip (Seq.length segment)
        Some(segment, rest)
    )

let chopStr segmentsize (source:string) =
    seq {
        for i in [0..segmentsize..source.Length-1] do
            yield source.Substring(i,segmentsize)
    }

let decoder =
    @"UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G".Split([|' ';'\n';'\r'|],StringSplitOptions.RemoveEmptyEntries)
    |> chop 2
    |> Seq.map (fun s->Array.ofSeq s |> fun x->x.[0],x.[1])
    |> Map.ofSeq
           
let decode (s:string) =
    chopStr 3 s
    |> Seq.map (fun x -> x|> Seq.map string|>String.concat "")    
    |> Seq.map (fun x -> decoder.[x])
    |> Seq.takeWhile (fun x -> x<>"Stop")
    |> String.concat ""
            
getData "prot" |> decode|>printfn "%s" 

// 10. SUBS. Finding a Motif in DNA
let rec allOccurrences (s:String) (t:String) =
    let ix=s.IndexOf(t)
    if (ix>=0) then 
        (ix+1)::(allOccurrences (s.Substring(ix+1)) t)
    else
        []
    |> Seq.scan (fun s t -> s+t) 0
    |> Seq.skip 1
    |> List.ofSeq   


[<EntryPoint>]
let main argv = 
    File.ReadAllText(argv.[0])
    |> fun x -> x.Trim()
    |> Seq.countBy (fun x->x)
    |> Map.ofSeq
    |> fun m -> printfn "%d %d %d %d" (guard m 'A') (guard m 'C') (guard m 'G') (guard m 'T')
    0 // return an integer exit code





// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.
open System
open System.Collections.Generic

// nwck
type Nwck =
    | OPEN
    | CLOSE
    | COMMA
    | Label of string
    | END

let rec lex (str:string) =
    if str="" then []
    else 
    match str.[0] with 
    | '(' -> OPEN :: (lex (str.Substring(1)))
    | ')' -> CLOSE :: (lex (str.Substring(1)))
    | ';' -> END :: (lex (str.Substring(1)))
    | ',' -> COMMA :: (lex (str.Substring(1)))
    | _ -> 
        let c = str.[0]
        let rest=lex (str.Substring(1))
        match rest with
        | Label(s) :: tail -> Label((string c) + s) :: tail
        | head :: tail -> Label(string c) :: head :: tail
        | [] -> Label(string c) :: []

type TreeState=(int*Dictionary<string,string list>)
let parse t (n:Nwck list) =
    match n with
    | [] -> t
    | head :: tail ->
        match head with
        | END -> t
        | 

// for now a tree is just an adjacency list
let rec parseNwck t str =
    match str.[0] with 
    | '(' ->
        let (l,label, remainder) = parseNodeList str 

    

[<EntryPoint>]
let main argv = 
    printfn "%A" argv
    0 // return an integer exit code

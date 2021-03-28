export interface Handler {
    name_token: string;
    handle: (parent: HTMLElement, body_raw: string) => void;
}
export interface InlineHandler {
    pattern: RegExp;
    make_span: (body_raw: string) => HTMLElement;
}

import { Handler, InlineHandler } from './handler';
export declare class MathHandler implements Handler {
    name_token: string;
    index_global: number;
    handle: (parent: HTMLElement, body_raw: string) => void;
}
export declare class InlineMathHandler implements InlineHandler {
    pattern: RegExp;
    make_span: (body_raw: string) => HTMLSpanElement;
}

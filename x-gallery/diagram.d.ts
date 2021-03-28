import { Handler } from './handler';
export declare const reset: () => void;
export declare class ImageDiagramHandler implements Handler {
    name_token: string;
    handle: (parent: HTMLElement, body_raw: string) => void;
}
export declare class P5DiagramHandler implements Handler {
    name_token: string;
    posts: any[];
    handle: (parent: HTMLElement, body_raw: string) => void;
}
export declare class ShaderDiagramHandler implements Handler {
    name_token: string;
    posts: any[];
    handle: (parent: HTMLElement, body_raw: string) => void;
}

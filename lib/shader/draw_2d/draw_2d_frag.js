window.shaders = window.shaders || {};
window.shaders.FragmentShader = `
#ifdef GL_ES
precision mediump float;//将精度设置为中等精度
#endif

varying vec2 v_position;

vec3 draw_circle(float margin, float radius, vec2 st){
    float pct = smoothstep(radius+margin, radius, length(st-vec2(0.5)))-smoothstep(radius, radius-margin, length(st-vec2(0.5)));
    return mix(vec3(1.0, 1.0, 1.0), vec3(0.4824, 0.3843, 0.7529), pct);
}

void main() {

    vec2 st = v_position * 0.5 + 0.5;

    vec3 color = vec3(0);

    st = fract(5.0*st);

    // Plot a line
    float margin = 0.01;
    float radius = 0.35;
    color = draw_circle(margin, radius, st);

    gl_FragColor = vec4(color, 1.0);
}
`;